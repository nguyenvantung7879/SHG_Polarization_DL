clc;
clear;
close all;
format long e;
warning off;
addpath("function\")

baseFolder = 'data_test\';
polarizerBase = fullfile(baseFolder, 'Polarizer');
angleBase     = fullfile(baseFolder, 'Angle');

%% =================== OUTPUT SETTINGS ===================
% Keep all porcine/tissue outputs in one folder, similar to tendon evaluation.
outRoot = fullfile(baseFolder, 'Tissue_Model_Evaluation_Results');
if ~isfolder(outRoot)
    mkdir(outRoot);
end

% Save predicted 2-page TIFF angle maps.
% Page 1 = predicted fiber orientation angle
% Page 2 = predicted peptide-pitch angle
writePredictionTifs = true;
predRoot = fullfile(outRoot, 'Prediction_TIFs');
if writePredictionTifs && ~isfolder(predRoot)
    mkdir(predRoot);
end

% Save representative Fig. 6 plot as PNG and FIG.
saveRepresentativeFigure = true;
figOutDir = fullfile(outRoot, 'Representative_Plots');
if saveRepresentativeFigure && ~isfolder(figOutDir)
    mkdir(figOutDir);
end

%% =================== Select sample order ===================
% Change this order to control row order in Fig. 6.
% Folder names must match the subfolders under Polarizer and Angle.
sampleOrder = {'Col I', 'Col II', 'Regerative tissue'};

% Select which file to show from each folder.
% Option A: use file index after alphabetical sorting.
fileIndexBySample = [10, 10, 1];

% Option B: use keyword matching. Leave '' to use fileIndexBySample.
% Example: fileKeywordBySample = {'ROI01', 'ROI03', 'week3'};
fileKeywordBySample = {'', '', ''};

%% =================== Load Model ====================
load('model\Finetune_Proposed_64RealPatch-Epoch-200-Selected32Patch-64-TrainPatchSize-16-Balance-ColI50-ColII50-Seed-7-Time-260.mat', 'net');
Net_Finetune = net;
clear net;

%% =================== Parameters =====================
% Common threshold/settings used for ALL rows.
% No sample-specific threshold/display cleanup here, so Col I / Col II / Reg
% are processed with the same settings for checking.
maskThreshold = 100;

% Display-only cleanup based on summed intensity percentile.
% 0 = disabled for all samples.
% Try 10, 15, 25 if you want to remove weak-signal pixels for all rows.
displayLowPercentile = 0;

% Display intensity normalization percentile.
% Robust display normalization for intensity.
% Lower the upper percentile if a few saturated white pixels make the image too dark.
% Recommended: [0.5, 98.5] or [0.5, 99].
intensityNormPercentile = [0.5, 98.5];

% Gamma < 1 brightens mid/low intensities after percentile normalization.
intensityGamma = 0.65;

% Histogram title / peak label display.
% These are applied directly inside plot_histogram_peak_only.
peakTitleFontSize = 14;
peakTitleFontWeight = 'bold';
peakLineWidth = 1.5;

% Export descriptive tables.
% IMPORTANT:
%   Selected figure peaks are captured directly from the plotted histograms,
%   so the CSV selected rows must match the displayed figure.
exportDescriptiveTables = true;
perFileOutFile         = fullfile(outRoot, 'Fig6_ColI_ColII_Reg_per_file.csv');
summaryBySampleOutFile = fullfile(outRoot, 'Fig6_ColI_ColII_Reg_summary_by_sample.csv');
figurePeakOutFile      = fullfile(outRoot, 'Fig6_selected_figure_peaks.csv');

% Display parameters from your old all-in-one script.
low = 1;
high = 99;
scale = 0.3;

% Match the old visual style. If the color direction is reversed, change to:
%   controlColors = jet(256);
controlColors = flipud(jet(256));

%% =================== Load selected files + predict ===================
SHG_all = cell(1, numel(sampleOrder));
Fiber_all = cell(1, numel(sampleOrder));
Pitch_all = cell(1, numel(sampleOrder));
selectedFiles = strings(1, numel(sampleOrder));

for s = 1:numel(sampleOrder)
    sampleName = sampleOrder{s};

    pFolder = fullfile(polarizerBase, sampleName);
    aFolder = fullfile(angleBase, sampleName);

    if ~isfolder(pFolder)
        error('Polarizer folder does not exist: %s', pFolder);
    end
    if ~isfolder(aFolder)
        warning('Angle folder does not exist, but prediction plot can still run: %s', aFolder);
    end

    [pshgPath, selectedFile] = selectTiffFile(pFolder, fileIndexBySample(s), fileKeywordBySample{s});
    selectedFiles(s) = string(selectedFile);

    fprintf('[%d/%d] %s | %s\n', s, numel(sampleOrder), sampleName, selectedFile);

    volume = tiffreadVolume(pshgPath);
    volume = single(volume);

    % Polarizer TIFF should contain only the P-SHG stack.
    % If it still contains extra pages, keep the first 18 polarization states.
    if size(volume, 3) >= 18
        SHG = volume(:, :, 1:18);
    else
        SHG = volume;
        warning('File %s has only %d pages. Using all pages as P-SHG stack.', selectedFile, size(volume,3));
    end

  
    currentThreshold = maskThreshold;

    sumIntensity = sum(SHG, 3);
    mask = sumIntensity <= currentThreshold;
    SHG(mask) = 0;

    [FiberPred, PitchPred, ~] = predict_image(Net_Finetune, SHG, currentThreshold);
    FiberPred(mask) = 0;
    PitchPred(mask) = 0;

    % Save selected predicted angle maps as 2-page TIFFs.
    if writePredictionTifs
        safeSampleName = makeSafeFolderName(sampleName);
        predSampleDir = fullfile(predRoot, safeSampleName);
        if ~isfolder(predSampleDir)
            mkdir(predSampleDir);
        end
        writeTwoPageAngleTiff(fullfile(predSampleDir, selectedFile), FiberPred, PitchPred);
    end

    SHG_all{s} = SHG;
    Fiber_all{s} = FiberPred;
    Pitch_all{s} = PitchPred;
end

fprintf('\nCommon settings:\n');
fprintf('  maskThreshold = %.3f\n', maskThreshold);
fprintf('  displayLowPercentile = %.3f\n', displayLowPercentile);
fprintf('  intensityNormPercentile = [%.3f %.3f]\n', intensityNormPercentile(1), intensityNormPercentile(2));
fprintf('  intensityGamma = %.3f\n', intensityGamma);

fprintf('\nSelected files:\n');
for s = 1:numel(sampleOrder)
    fprintf('  %s: %s\n', sampleOrder{s}, selectedFiles(s));
end
fprintf('\n');

%% =================== Fig. 6 visualization ===================
[figureFiberPeaks, figurePitchPeaks] = plot_fig6_all_in_one_results( ...
    SHG_all{1}, Fiber_all{1}, Pitch_all{1}, ...
    SHG_all{2}, Fiber_all{2}, Pitch_all{2}, ...
    SHG_all{3}, Fiber_all{3}, Pitch_all{3}, ...
    low, high, scale, controlColors, ...
    displayLowPercentile, intensityNormPercentile, intensityGamma, ...
    peakTitleFontSize, peakTitleFontWeight, peakLineWidth);

if saveRepresentativeFigure
    figHandles = findall(0, 'Type', 'figure');
    figHandles = flipud(figHandles(:));

    for ff = 1:numel(figHandles)
        pngPath = fullfile(figOutDir, sprintf('Fig6_all_in_one_selected_order_%02d.png', ff));
        figPath = fullfile(figOutDir, sprintf('Fig6_all_in_one_selected_order_%02d.fig', ff));

        try
            saveas(figHandles(ff), pngPath);
            saveas(figHandles(ff), figPath);
        catch ME
            warning('Could not save representative figure %d: %s', ff, ME.message);
        end
    end

    fprintf('Saved representative plots to: %s\n', figOutDir);
end

FigurePeakTable = table(string(sampleOrder(:)), selectedFiles(:), ...
    double(figureFiberPeaks(:)), double(figurePitchPeaks(:)), ...
    'VariableNames', {'Sample','SelectedFile','FiberPredPeak_deg','PitchPredPeak_deg'});

if exportDescriptiveTables
    writetable(FigurePeakTable, figurePeakOutFile);
    fprintf('Saved selected figure peaks: %s\n', figurePeakOutFile);

    [PerFileTable, SummaryBySample] = build_fig6_all_file_descriptive_tables( ...
        sampleOrder, selectedFiles, FigurePeakTable, ...
        polarizerBase, angleBase, ...
        Net_Finetune, maskThreshold, displayLowPercentile, ...
        writePredictionTifs, predRoot);

    writetable(PerFileTable, perFileOutFile);
    writetable(SummaryBySample, summaryBySampleOutFile);

    fprintf('Saved per-file table: %s\n', perFileOutFile);
    fprintf('Saved summary-by-sample table: %s\n\n', summaryBySampleOutFile);
end

% Optional save
% saveas(gcf, 'Fig6_all_in_one_selected_order.png');
% saveas(gcf, 'Fig6_all_in_one_selected_order.fig');


%% =================== Descriptive table helper functions ===================
function [PerFileTable, SummaryBySample] = build_fig6_all_file_descriptive_tables( ...
    sampleOrder, selectedFiles, FigurePeakTable, polarizerBase, angleBase, ...
    Net_Finetune, maskThreshold, displayLowPercentile, writePredictionTifs, predRoot)

    PerFileTable = table();

    for s = 1:numel(sampleOrder)
        sampleName = sampleOrder{s};
        pFolder = fullfile(polarizerBase, sampleName);
        aFolder = fullfile(angleBase, sampleName);

        fileList = listTiffFilesSorted(pFolder);
        fprintf('Summary scan | %s: %d TIFF files\n', sampleName, numel(fileList));

        for f = 1:numel(fileList)
            fileName = fileList(f).name;
            pshgPath = fullfile(pFolder, fileName);

            isFigureSelected = strcmp(string(fileName), string(selectedFiles(s)));

            try
                row = evaluate_one_fig6_file_for_table( ...
                    sampleName, f, fileName, pshgPath, aFolder, Net_Finetune, ...
                    maskThreshold, displayLowPercentile, isFigureSelected, ...
                    writePredictionTifs, predRoot);

                % Force selected rows to use the exact peaks shown in Fig. 6.
                % This is the most reliable way to guarantee CSV/figure matching.
                if isFigureSelected
                    row.FiberPredPeak_deg = FigurePeakTable.FiberPredPeak_deg(s);
                    row.PitchPredPeak_deg = FigurePeakTable.PitchPredPeak_deg(s);

                    if isfinite(row.FiberTheoryPeak_deg)
                        row.FiberPeakAbsError_deg = circularAbsDiff180(row.FiberPredPeak_deg, row.FiberTheoryPeak_deg);
                    end
                    if isfinite(row.PitchTheoryPeak_deg)
                        row.PitchPeakAbsError_deg = abs(row.PitchPredPeak_deg - row.PitchTheoryPeak_deg);
                    end
                end

                PerFileTable = [PerFileTable; row]; %#ok<AGROW>
            catch ME
                warning('Failed descriptive table for %s | %s: %s', sampleName, fileName, ME.message);
            end
        end
    end

    SummaryBySample = summarize_fig6_per_file_table(PerFileTable);
end

function row = evaluate_one_fig6_file_for_table( ...
    sampleName, fileIndex, fileName, pshgPath, aFolder, Net_Finetune, ...
    maskThreshold, displayLowPercentile, isFigureSelected, writePredictionTifs, predRoot)

    volume = tiffreadVolume(pshgPath);
    volume = single(volume);

    if size(volume, 3) >= 18
        SHG = volume(:, :, 1:18);
    else
        SHG = volume;
    end

    sumIntensityRaw = sum(SHG, 3);
    mask = sumIntensityRaw <= maskThreshold;
    SHG(mask) = 0;

    [FiberPred, PitchPred, ~] = predict_image(Net_Finetune, SHG, maskThreshold);

    % Save predicted angle maps for every tissue file scanned in the descriptive table.
    if writePredictionTifs
        safeSampleName = makeSafeFolderName(sampleName);
        predSampleDir = fullfile(predRoot, safeSampleName);
        if ~isfolder(predSampleDir)
            mkdir(predSampleDir);
        end
        writeTwoPageAngleTiff(fullfile(predSampleDir, fileName), FiberPred, PitchPred);
    end

    [FiberPredClean, PitchPredClean, finalMask] = clean_maps_like_fig6_histogram( ...
        SHG, FiberPred, PitchPred, displayLowPercentile);

    sumIntensity = sum(SHG, 3);
    positiveVals = double(sumIntensity(sumIntensity > 0 & isfinite(sumIntensity)));
    if isempty(positiveVals)
        meanPSHG = NaN;
    else
        meanPSHG = mean(positiveVals, 'omitnan');
    end

    validMask = isfinite(FiberPredClean) & isfinite(PitchPredClean);
    validPixelFraction = nnz(validMask) / numel(validMask);

    fiberPredPeak = peakLikeFigureSafe(FiberPredClean, 1:180);
    pitchPredPeak = peakLikeFigureSafe(PitchPredClean, 40:60);

    fiberPredVals = clean_vec_for_hist(FiberPredClean, [1, 180]);
    pitchPredVals = clean_vec_for_hist(PitchPredClean, [40, 60]);

    fiberPredMean = circularMean180(fiberPredVals);
    fiberPredSD   = circularStd180(fiberPredVals);
    pitchPredMean = mean(pitchPredVals, 'omitnan');
    pitchPredSD   = std(pitchPredVals, 'omitnan');

    [anglePath, hasTheory] = findMatchingAngleTiff(aFolder, fileName);

    fiberTheoryPeak = NaN;
    pitchTheoryPeak = NaN;
    fiberPeakAbsError = NaN;
    pitchPeakAbsError = NaN;

    if hasTheory
        infoA = imfinfo(anglePath);
        if numel(infoA) >= 2
            FiberTheory = double(imread(anglePath, 'Index', 1));
            PitchTheory = double(imread(anglePath, 'Index', 2));

            FiberTheory(finalMask) = NaN;
            PitchTheory(finalMask) = NaN;
            FiberTheory(FiberTheory <= 0 | FiberTheory >= 180) = NaN;
            PitchTheory(PitchTheory <= 0) = NaN;

            fiberTheoryPeak = peakLikeFigureSafe(FiberTheory, 1:180);
            pitchTheoryPeak = peakLikeFigureSafe(PitchTheory, 40:60);

            fiberPeakAbsError = circularAbsDiff180(fiberPredPeak, fiberTheoryPeak);
            pitchPeakAbsError = abs(pitchPredPeak - pitchTheoryPeak);
        end
    end

    if contains(lower(string(sampleName)), 'reg')
        noteText = "Heterogeneous biological sample; descriptive only";
    else
        noteText = "Representative collagen sample; descriptive only";
    end

    row = table( ...
        string(sampleName), double(fileIndex), string(fileName), logical(isFigureSelected), ...
        double(meanPSHG), double(validPixelFraction), ...
        double(fiberPredPeak), double(fiberPredMean), double(fiberPredSD), ...
        double(pitchPredPeak), double(pitchPredMean), double(pitchPredSD), ...
        double(fiberTheoryPeak), double(pitchTheoryPeak), ...
        double(fiberPeakAbsError), double(pitchPeakAbsError), string(noteText), ...
        'VariableNames', {'Sample','FileIndex','SelectedFile','IsFigureSelected', ...
        'MeanPSHGIntensity','ValidPixelFraction', ...
        'FiberPredPeak_deg','FiberPredMean_deg','FiberPredSD_deg', ...
        'PitchPredPeak_deg','PitchPredMean_deg','PitchPredSD_deg', ...
        'FiberTheoryPeak_deg','PitchTheoryPeak_deg', ...
        'FiberPeakAbsError_deg','PitchPeakAbsError_deg','Note'});
end

function [FiberClean, PitchClean, finalMask] = clean_maps_like_fig6_histogram(SHG, FiberPred, PitchPred, displayLowPercentile)
    FiberClean = double(gather(FiberPred));
    PitchClean = double(gather(PitchPred));

    baseMask = compute_hist_mask_from_pshg(SHG);
    weakMask = compute_low_intensity_percentile_mask(SHG, displayLowPercentile);
    finalMask = baseMask | weakMask;

    FiberClean(finalMask) = NaN;
    PitchClean(finalMask) = NaN;

    FiberClean(FiberClean <= 0 | FiberClean >= 180) = NaN;
    PitchClean(PitchClean <= 0) = NaN;
end

function peakVal = peakLikeFigureSafe(img, rangeVals)
    x = clean_vec_for_hist(img, [min(rangeVals), max(rangeVals)]);
    if isempty(x)
        peakVal = NaN;
    else
        peakVal = safeHistogramPeakForPlot(x, rangeVals);
    end
end

function peakVal = safeHistogramPeakForPlot(x, rangeVals)
    % Try Gaussian peak first. If invalid/negative/out-of-range, fallback to histogram mode.
    x = double(x(:));
    x = x(isfinite(x));

    minR = min(rangeVals);
    maxR = max(rangeVals);
    x = x(x >= minR & x <= maxR);

    if isempty(x)
        peakVal = NaN;
        return;
    end

    peakVal = NaN;

    try
        [~, ~, fitPeak, ~] = fitGaussianHistogram(x, rangeVals, 'angle', 0, []);
        if isfinite(fitPeak) && fitPeak >= minR && fitPeak <= maxR
            peakVal = fitPeak;
        end
    catch
        peakVal = NaN;
    end

    if ~isfinite(peakVal)
        peakVal = histogramModePeak(x, rangeVals);
    end
end

function peakVal = histogramModePeak(x, rangeVals)
    x = double(x(:));
    x = x(isfinite(x));

    if isempty(x)
        peakVal = NaN;
        return;
    end

    rangeVals = double(rangeVals(:)).';
    if numel(rangeVals) < 2
        peakVal = NaN;
        return;
    end

    binStep = median(diff(rangeVals));
    edges = rangeVals;

    if max(x) >= max(edges)
        edges = [edges, max(edges) + binStep];
    end

    counts = histcounts(x, edges);
    if isempty(counts) || max(counts) <= 0
        peakVal = NaN;
        return;
    end

    centers = (edges(1:end-1) + edges(2:end)) ./ 2;
    [~, idx] = max(counts);
    peakVal = centers(idx);
end

function SummaryBySample = summarize_fig6_per_file_table(T)
    sampleList = unique(T.Sample, 'stable');
    SummaryBySample = table();

    for i = 1:numel(sampleList)
        sampleName = sampleList(i);
        idx = T.Sample == sampleName;

        row = table();
        row.Sample = sampleName;
        row.N = sum(idx);

        row.MeanPSHGIntensity = mean(T.MeanPSHGIntensity(idx), 'omitnan');
        row.MeanPSHGIntensity_SD = std(T.MeanPSHGIntensity(idx), 'omitnan');

        row.ValidPixelFraction = mean(T.ValidPixelFraction(idx), 'omitnan');
        row.ValidPixelFraction_SD = std(T.ValidPixelFraction(idx), 'omitnan');

        row.FiberPredPeak_deg = circularMean180(T.FiberPredPeak_deg(idx));
        row.FiberPredPeak_deg_SD = circularStd180(T.FiberPredPeak_deg(idx));

        row.PitchPredPeak_deg = mean(T.PitchPredPeak_deg(idx), 'omitnan');
        row.PitchPredPeak_deg_SD = std(T.PitchPredPeak_deg(idx), 'omitnan');

        row.FiberTheoryPeak_deg = circularMean180(T.FiberTheoryPeak_deg(idx));
        row.FiberTheoryPeak_deg_SD = circularStd180(T.FiberTheoryPeak_deg(idx));

        row.PitchTheoryPeak_deg = mean(T.PitchTheoryPeak_deg(idx), 'omitnan');
        row.PitchTheoryPeak_deg_SD = std(T.PitchTheoryPeak_deg(idx), 'omitnan');

        row.FiberPeakAbsError_deg = mean(T.FiberPeakAbsError_deg(idx), 'omitnan');
        row.FiberPeakAbsError_deg_SD = std(T.FiberPeakAbsError_deg(idx), 'omitnan');

        row.PitchPeakAbsError_deg = mean(T.PitchPeakAbsError_deg(idx), 'omitnan');
        row.PitchPeakAbsError_deg_SD = std(T.PitchPeakAbsError_deg(idx), 'omitnan');

        row.Note = T.Note(find(idx, 1, 'first'));

        SummaryBySample = [SummaryBySample; row]; %#ok<AGROW>
    end
end

function fileList = listTiffFilesSorted(folderPath)
    fileList = [dir(fullfile(folderPath, '*.tif')); dir(fullfile(folderPath, '*.tiff'))];

    if isempty(fileList)
        error('No TIFF files found in: %s', folderPath);
    end

    [~, sortIdx] = sort({fileList.name});
    fileList = fileList(sortIdx);
end

function [anglePath, hasTheory] = findMatchingAngleTiff(angleFolder, selectedFile)
    hasTheory = false;
    anglePath = '';

    if ~isfolder(angleFolder)
        return;
    end

    [~, sampleName, ext] = fileparts(selectedFile);

    candidate1 = fullfile(angleFolder, selectedFile);
    if isfile(candidate1)
        anglePath = candidate1;
        hasTheory = true;
        return;
    end

    if strcmpi(ext, '.tif')
        candidate2 = fullfile(angleFolder, [sampleName, '.tiff']);
    else
        candidate2 = fullfile(angleFolder, [sampleName, '.tif']);
    end

    if isfile(candidate2)
        anglePath = candidate2;
        hasTheory = true;
    end
end

function d = circularAbsDiff180(a, b)
    if ~isfinite(a) || ~isfinite(b)
        d = NaN;
        return;
    end
    d = abs(mod((a - b) + 90, 180) - 90);
end

function mu = circularMean180(x)
    x = double(x(:));
    x = x(isfinite(x) & x > 0 & x < 180);

    if isempty(x)
        mu = NaN;
        return;
    end

    theta = deg2rad(2*x);
    c = mean(cos(theta));
    s = mean(sin(theta));
    mu = mod(rad2deg(atan2(s, c)) / 2, 180);
end

function sd = circularStd180(x)
    x = double(x(:));
    x = x(isfinite(x) & x > 0 & x < 180);

    if isempty(x)
        sd = NaN;
        return;
    end

    theta = deg2rad(2*x);
    R = sqrt(mean(cos(theta))^2 + mean(sin(theta))^2);
    R = max(min(R, 1), eps);
    sd = rad2deg(sqrt(-2*log(R))) / 2;
end


%% =================== File selection helper ===================
function [filePath, fileName] = selectTiffFile(folderPath, fileIndex, keyword)
    fileList = [dir(fullfile(folderPath, '*.tif')); dir(fullfile(folderPath, '*.tiff'))];

    if isempty(fileList)
        error('No TIFF files found in: %s', folderPath);
    end

    [~, sortIdx] = sort({fileList.name});
    fileList = fileList(sortIdx);

    if nargin >= 3 && ~isempty(keyword)
        hit = false(numel(fileList), 1);
        for i = 1:numel(fileList)
            hit(i) = contains(lower(fileList(i).name), lower(keyword));
        end
        fileList = fileList(hit);

        if isempty(fileList)
            error('No TIFF file containing keyword "%s" in: %s', keyword, folderPath);
        end

        fileIndex = 1;
    end

    if fileIndex < 1 || fileIndex > numel(fileList)
        error('fileIndex %d is out of range. Folder %s has %d TIFF files.', fileIndex, folderPath, numel(fileList));
    end

    fileName = fileList(fileIndex).name;
    filePath = fullfile(folderPath, fileName);
end



function writeTwoPageAngleTiff(outPath, fiberMap, pitchMap)
    % Write a 2-page 32-bit float TIFF:
    %   page 1 = fiber orientation angle
    %   page 2 = peptide-pitch angle
    fiberMap = single(gather(fiberMap));
    pitchMap = single(gather(pitchMap));

    if isfile(outPath)
        delete(outPath);
    end

    t = Tiff(outPath, 'w');

    tagstruct = makeFloatTiffTag(fiberMap);
    setTag(t, tagstruct);
    write(t, fiberMap);

    writeDirectory(t);
    tagstruct = makeFloatTiffTag(pitchMap);
    setTag(t, tagstruct);
    write(t, pitchMap);

    close(t);
end

function tagstruct = makeFloatTiffTag(img)
    tagstruct.ImageLength = size(img, 1);
    tagstruct.ImageWidth = size(img, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.Software = 'MATLAB';
end

function safeName = makeSafeFolderName(nameIn)
    safeName = char(string(nameIn));
    safeName = strrep(safeName, ' ', '_');
    safeName = strrep(safeName, '/', '_');
    safeName = strrep(safeName, '\', '_');
end

%% =================== Local plotting function ===================
function [figureFiberPeaks, figurePitchPeaks] = plot_fig6_all_in_one_results( ...
    SHG_1, Fiber_1, Pitch_1, ...
    SHG_2, Fiber_2, Pitch_2, ...
    SHG_3, Fiber_3, Pitch_3, ...
    low, high, scale, controlColors, ...
    displayLowPercentile, intensityNormPercentile, intensityGamma, ...
    peakTitleFontSize, peakTitleFontWeight, peakLineWidth)

    % -------- clean data --------
    SHG_all   = {double(gather(SHG_1)),   double(gather(SHG_2)),   double(gather(SHG_3))};
    Fiber_all = {double(gather(Fiber_1)), double(gather(Fiber_2)), double(gather(Fiber_3))};
    Pitch_all = {double(gather(Pitch_1)), double(gather(Pitch_2)), double(gather(Pitch_3))};

    for k = 1:3
        mask = compute_hist_mask_from_pshg(SHG_all{k});

        % Display-only cleanup: remove weak summed-intensity pixels.
        % This is especially useful for regenerative tissue where weak SHG
        % pixels can produce scattered peptide-pitch colors.
        pLowMask = displayLowPercentile;

        weakMask = compute_low_intensity_percentile_mask(SHG_all{k}, pLowMask);
        mask = mask | weakMask;

        Fiber_all{k}(mask) = NaN;
        Pitch_all{k}(mask) = NaN;

        Fiber_all{k}(Fiber_all{k} <= 0 | Fiber_all{k} >= 180) = NaN;
        Pitch_all{k}(Pitch_all{k} <= 0) = NaN;
    end

    map_fiber = controlColors;
    map_pitch = controlColors;

    fiber_clim = [0 180];
    pitch_clim = [40 60];

    % Peaks actually displayed in the figure.
    figureFiberPeaks = NaN(3, 1);
    figurePitchPeaks = NaN(3, 1);

    % -------- overlays: use coat_SHG_with_angle if available --------
    fiber_rgb = cell(1,3);
    for k = 1:3
        fiber_rgb{k} = make_fiber_overlay_rgb(SHG_all{k}, Fiber_all{k}, map_fiber, low, high, scale);
    end

    % -------- figure layout: white background, 3 rows x 5 columns --------
    fig = figure('Units', 'pixels', ...
                 'Position', [10, 40, 2048, 950], ...
                 'Color', 'w', ...
                 'Name', 'Fig. 6 all-in-one old style white background');

    row_h = 0.245;
    row_gap = 0.060;
    y0 = 0.055;
    y_pos = y0 + (2:-1:0) .* (row_h + row_gap);

    % Columns. Layout is tuned to mimic your old result:
    % intensity | fiber overlay | fiber histogram | pitch map | pitch histogram
    x_int  = 0.015; w_int  = 0.130;
    x_fib  = 0.205; w_fib  = 0.130;
    x_fhis = 0.420; w_fhis = 0.150;
    x_pit  = 0.630; w_pit  = 0.130;
    x_phis = 0.825; w_phis = 0.150;

    % Full-height colorbars
    cb_w = 0.012;
    cb_int_x = 0.153;
    cb_fib_x = 0.343;
    cb_pit_x = 0.768;
    cb_y = y0;
    cb_h = 3*row_h + 2*row_gap;

    for r = 1:3
        y = y_pos(r);

        % Column 1: P-SHG intensity
        ax1 = axes('Parent', fig, 'Position', [x_int, y, w_int, row_h]);
        intensityImg = compute_pshg_intensity(SHG_all{r}, intensityNormPercentile, intensityGamma);
        imagesc(ax1, intensityImg, [0 1]);
        axis(ax1, 'image', 'off');
        colormap(ax1, gray(256));
        set(ax1, 'Color', 'k');
        draw_panel_box(ax1, intensityImg, 'w');

        % Column 2: fiber orientation overlay on SHG intensity
        ax2 = axes('Parent', fig, 'Position', [x_fib, y, w_fib, row_h]);
        plot_fiber_rgb_no_roi(ax2, fiber_rgb{r});

        % Column 3: fiber histogram, peak line only
        ax3 = axes('Parent', fig, 'Position', [x_fhis, y, w_fhis, row_h]);
        figureFiberPeaks(r) = plot_histogram_peak_only(ax3, Fiber_all{r}, 1:180, 'Fiber orientation (\circ)', ...
            peakTitleFontSize, peakTitleFontWeight, peakLineWidth);
        set(ax3, 'FontSize', 9, 'LineWidth', 1);

        % Column 4: peptide-pitch angle map
        ax4 = axes('Parent', fig, 'Position', [x_pit, y, w_pit, row_h]);
        plot_angle_map(ax4, Pitch_all{r}, map_pitch, pitch_clim);

        % Column 5: peptide-pitch histogram, peak line only
        ax5 = axes('Parent', fig, 'Position', [x_phis, y, w_phis, row_h]);
        figurePitchPeaks(r) = plot_histogram_peak_only(ax5, Pitch_all{r}, 40:60, 'Peptide-pitch angle (\circ)', ...
            peakTitleFontSize, peakTitleFontWeight, peakLineWidth);
        set(ax5, 'FontSize', 9, 'LineWidth', 1);
    end

    % -------- full-height colorbars --------
    cb_ax0 = axes('Parent', fig, 'Position', [cb_int_x, cb_y, cb_w, cb_h], 'Visible', 'off', 'Color', 'w');
    colormap(cb_ax0, gray(256));
    caxis(cb_ax0, [0 1]);
    cb0 = colorbar(cb_ax0, 'Location', 'eastoutside');
    cb0.Units = 'normalized';
    cb0.Position = [cb_int_x, cb_y, cb_w, cb_h];
    cb0.Ticks = 0:0.1:1;
    cb0.FontSize = 14;
    cb0.Color = 'k';
    cb0.Label.String = 'P-SHG intensity';
    cb0.Label.FontSize = 14;
    cb0.Label.Color = 'k';

    cb_ax1 = axes('Parent', fig, 'Position', [cb_fib_x, cb_y, cb_w, cb_h], 'Visible', 'off', 'Color', 'w');
    colormap(cb_ax1, map_fiber);
    caxis(cb_ax1, fiber_clim);
    cb1 = colorbar(cb_ax1, 'Location', 'eastoutside');
    cb1.Units = 'normalized';
    cb1.Position = [cb_fib_x, cb_y, cb_w, cb_h];
    cb1.Ticks = 0:30:180;
    cb1.FontSize = 14;
    cb1.Color = 'k';
    cb1.Label.String = 'Fiber orientation (\circ)';
    cb1.Label.FontSize = 14;
    cb1.Label.Color = 'k';

    cb_ax2 = axes('Parent', fig, 'Position', [cb_pit_x, cb_y, cb_w, cb_h], 'Visible', 'off', 'Color', 'w');
    colormap(cb_ax2, map_pitch);
    caxis(cb_ax2, pitch_clim);
    cb2 = colorbar(cb_ax2, 'Location', 'eastoutside');
    cb2.Units = 'normalized';
    cb2.Position = [cb_pit_x, cb_y, cb_w, cb_h];
    cb2.Ticks = 40:5:60;
    cb2.FontSize = 14;
    cb2.Color = 'k';
    cb2.Label.String = 'Peptide-pitch angle (\circ)';
    cb2.Label.FontSize = 14;
    cb2.Label.Color = 'k';
end

%% =================== Helper functions ===================
function mask = compute_hist_mask_from_pshg(I_SHG)
    % Mask based on summed P-SHG intensity across polarization states.
    I_SHG = double(gather(I_SHG));

    if ndims(I_SHG) == 2
        sumIntensity = I_SHG;
    else
        sumIntensity = sum(I_SHG, 3);
    end

    mask = sumIntensity <= 0;
    mask = logical(mask);
end

function mask = compute_low_intensity_percentile_mask(I_SHG, pLow)
    % Display-only mask based on summed P-SHG intensity percentile.
    % pLow = 0 disables this extra cleanup.
    I_SHG = double(gather(I_SHG));

    if nargin < 2 || isempty(pLow) || pLow <= 0
        if ndims(I_SHG) == 2
            mask = false(size(I_SHG));
        else
            mask = false(size(I_SHG,1), size(I_SHG,2));
        end
        return;
    end

    if ndims(I_SHG) == 2
        sumIntensity = I_SHG;
    else
        sumIntensity = sum(I_SHG, 3);
    end

    positiveVals = sumIntensity(sumIntensity > 0 & isfinite(sumIntensity));
    if isempty(positiveVals)
        mask = true(size(sumIntensity));
        return;
    end

    cutoff = prctile(positiveVals, pLow);
    mask = sumIntensity <= cutoff;
end

function rgb = make_fiber_overlay_rgb(I_SHG, angleMap, cmap, low, high, scale)
    % Prefer your original coating function when it exists in MATLAB path.
    if exist('coat_SHG_with_angle', 'file') == 2
        try
            tmp = coat_SHG_with_angle(I_SHG, angleMap, cmap, low, high, scale);
            rgb = im2double(tmp(:,:,:,min(9, size(tmp,4))));
            return;
        catch
            % Fall back to local overlay if the external function fails.
        end
    end

    % Fallback overlay, only used when coat_SHG_with_angle is unavailable.
    I = compute_pshg_intensity(I_SHG, intensityNormPercentile_default(), 0.65);
    A = double(gather(angleMap));
    A = flipud(A);

    valid = isfinite(A) & A > 0 & A < 180;

    idx = round(A ./ 180 .* (size(cmap,1)-1)) + 1;
    idx(~valid) = 1;
    idx = max(min(idx, size(cmap,1)), 1);

    angleRGB = ind2rgb(idx, cmap);

    I2 = I;
    positiveVals = I2(I2 > 0);
    if ~isempty(positiveVals)
        lo = prctile(positiveVals, max(0, min(100, low)));
        hi = prctile(positiveVals, max(0, min(100, high)));
    else
        lo = min(I2(:)); hi = max(I2(:)) + eps;
    end

    if hi <= lo
        lo = min(I2(:)); hi = max(I2(:)) + eps;
    end

    I2 = (I2 - lo) ./ (hi - lo + eps);
    I2 = max(min(I2, 1), 0);
    grayRGB = repmat(I2, [1 1 3]);

    rgb = (1-scale) .* grayRGB + scale .* angleRGB;
    rgb(~repmat(valid, [1 1 3])) = 0;
end

function peakVal = plot_histogram_peak_only(ax, img, rangeVals, xLabelStr, peakTitleFontSize, peakTitleFontWeight, peakLineWidth)
    % Histogram with red peak line only.
    % This function returns the exact peak displayed in the figure.

    if nargin < 5 || isempty(peakTitleFontSize)
        peakTitleFontSize = 12;
    end
    if nargin < 6 || isempty(peakTitleFontWeight)
        peakTitleFontWeight = 'bold';
    end
    if nargin < 7 || isempty(peakLineWidth)
        peakLineWidth = 1.5;
    end

    axes(ax); %#ok<LAXES>
    x = clean_vec_for_hist(img, [min(rangeVals), max(rangeVals)]);

    if isempty(x)
        histogram(ax, NaN, rangeVals, 'EdgeColor', 'none', 'FaceAlpha', 1);
        peakVal = NaN;
    else
        histogram(ax, x, rangeVals, 'EdgeColor', 'none', 'FaceAlpha', 1);
        peakVal = safeHistogramPeakForPlot(x, rangeVals);
    end

    hold(ax, 'on');

    if isfinite(peakVal)
        xline(ax, peakVal, 'r', 'LineWidth', peakLineWidth);
        peakText = ['Peak = ' num2str(peakVal, '%.2f') '\circ'];
    else
        peakText = 'Peak = NaN';
    end

    title(ax, '');

    text(ax, 0.5, 1.00, peakText, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', peakTitleFontSize, ...
        'FontWeight', peakTitleFontWeight, ...
        'Color', 'k', ...
        'Interpreter', 'tex', ...
        'Clipping', 'off');

    xlabel(ax, xLabelStr, 'FontSize', 9);
    ylabel(ax, 'Counts', 'FontSize', 9);
    xlim(ax, [min(rangeVals), max(rangeVals)]);
    axis(ax, 'square');
    box(ax, 'on');
    hold(ax, 'off');
end

function [fitresult, FWHM, Center, MaxCount] = fitGaussianHistogram(data, range, x_label, plot_figure, yMax)
    if nargin < 2 || isempty(range)
        dataTmp = double(data(:));
        dataTmp = dataTmp(isfinite(dataTmp));

        if isempty(dataTmp)
            fitresult = [];
            FWHM = NaN;
            Center = NaN;
            MaxCount = NaN;
            return;
        end

        range = floor(min(dataTmp)):1:ceil(max(dataTmp));
    end

    if nargin < 3 || isempty(x_label), x_label = 'angle'; end %#ok<NASGU>
    if nargin < 4 || isempty(plot_figure), plot_figure = 0; end
    if nargin < 5, yMax = []; end %#ok<NASGU>

    data = double(data(:));
    range = double(range(:)).';
    data = data(isfinite(data));

    if isempty(data) || numel(range) < 2
        fitresult = [];
        FWHM = NaN;
        Center = NaN;
        MaxCount = NaN;
        return;
    end

    X = double((range(1:end-1) + range(2:end)) ./ 2);
    Num_Fiber = double(histcounts(data, range));
    MaxCount = double(max(Num_Fiber));

    if MaxCount <= 0 || numel(X) < 3
        fitresult = [];
        FWHM = NaN;
        Center = NaN;
        return;
    end

    try
        ft_G = fittype('a*exp(-(x-b).^2 / c.^2 / 2)', 'independent', 'x', 'dependent', 'y');
        opts_G = fitoptions('Method', 'NonlinearLeastSquares');
        opts_G.Display = 'Off';
        opts_G.Lower = double([0 -inf 0]);

        [xData, yData] = prepareCurveData(double(X(:)), double(Num_Fiber(:)));
        xData = double(xData);
        yData = double(yData);

        [I_P, N_max] = max(Num_Fiber);
        opts_G.StartPoint = double([0.5 * I_P, X(N_max), 2]);

        [fitresult, ~] = fit(xData, yData, ft_G, opts_G);
        PC = double(coeffvalues(fitresult));
        FWHM = PC(3) * 2 * sqrt(2*log(2));
        Center = PC(2);
    catch
        fitresult = [];
        FWHM = NaN;
        [~, N_max] = max(Num_Fiber);
        Center = X(N_max);
    end

    if plot_figure == 1
        histogram(data, range, 'EdgeColor', 'none', 'FaceAlpha', 1); hold on;
        xline(Center, 'r', 'LineWidth', 1);
        title(['Peak = ' num2str(Center, '%.2f') '\circ']);
        xlabel('angle'); ylabel('Counts');
        xlim([min(range), max(range)]);
        axis square;
        hold off;
    end
end

function x = clean_vec_for_hist(img, climRange)
    x = double(gather(img(:)));
    x = x(isfinite(x));
    x = x(x >= climRange(1) & x <= climRange(2));
end

function plot_fiber_rgb_no_roi(ax, rgb_img)
    [nr, nc, ~] = size(rgb_img);
    imshow(rgb_img, 'Parent', ax);
    axis(ax, 'image', 'off');
    hold(ax, 'on');
    rectangle(ax, 'Position', [1, 1, nc-1, nr-1], ...
              'EdgeColor', 'w', 'LineWidth', 1.0, 'Clipping', 'on');
    hold(ax, 'off');
end

function plot_angle_map(ax, angleMap, cmap, climRange)
    angleMap = double(gather(angleMap));
    angleMap = flipud(angleMap);
    angleMap(angleMap <= climRange(1) | angleMap >= climRange(2)) = NaN;
    mask = isnan(angleMap);

    hold(ax, 'on');
    imagesc(ax, zeros([size(angleMap), 3]));
    h = imagesc(ax, angleMap);
    h.AlphaData = double(~mask);
    axis(ax, 'image', 'off');
    colormap(ax, cmap);
    caxis(ax, climRange);
    set(ax, 'Color', 'k', 'YDir', 'normal');
    rectangle(ax, 'Position', [0.5, 0.5, size(angleMap,2)-1, size(angleMap,1)-1], ...
              'EdgeColor', 'w', 'LineWidth', 1.0, 'Clipping', 'on');
    hold(ax, 'off');
end

function out = compute_pshg_intensity(Iin, normPercentile, gammaVal)
    % Display summed P-SHG intensity across polarization states.
    % Robust percentile clipping avoids a few saturated white pixels making
    % the whole image too dark.
    if nargin < 2 || isempty(normPercentile)
        normPercentile = [0.5, 98.5];
    end
    if nargin < 3 || isempty(gammaVal)
        gammaVal = 0.65;
    end

    Iin = double(gather(Iin));
    Iin(Iin < 0) = 0;

    if ndims(Iin) == 2
        out = Iin;
    else
        out = sum(Iin, 3);
        out = squeeze(out);
        while ndims(out) > 2
            out = sum(out, ndims(out));
        end
    end

    % Flip vertically so intensity matches fiber and pitch maps.
    % out = flipud(out);

    positiveVals = out(out > 0 & isfinite(out));
    if isempty(positiveVals)
        out = zeros(size(out));
        return;
    end

    % First remove extreme saturated outliers from the statistics only.
    satCut = prctile(positiveVals, 99.95);
    statVals = positiveVals(positiveVals <= satCut);
    if numel(statVals) < 20
        statVals = positiveVals;
    end

    lo = prctile(statVals, normPercentile(1));
    hi = prctile(statVals, normPercentile(2));

    if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
        lo = min(statVals);
        hi = max(statVals);
    end

    out = (out - lo) ./ (hi - lo + eps);
    out = max(min(out, 1), 0);

    % Brighten mid/low intensities for display.
    out = out .^ gammaVal;
end

function intensityNormPercentile = intensityNormPercentile_default()
    intensityNormPercentile = [0.5, 98.5];
end

function draw_panel_box(ax, img, colorSpec)
    hold(ax, 'on');
    rectangle(ax, 'Position', [0.5, 0.5, size(img,2)-1, size(img,1)-1], ...
              'EdgeColor', colorSpec, 'LineWidth', 1.0, 'Clipping', 'on');
    hold(ax, 'off');
end
