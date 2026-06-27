%% run_tendon_eval.m
clc;
clear;
close all;
format long e;
warning off;
addpath("function\")

%% ===================== USER SETTINGS =====================
baseFolder = 'Data_Test\';

pshgFolder  = fullfile(baseFolder, 'Polarizer', 'Tendon');
angleFolder = fullfile(baseFolder, 'Angle',     'Tendon');

outRoot = fullfile(baseFolder, 'Tendon_Model_Evaluation_Results');
if ~isfolder(outRoot), mkdir(outRoot); end

% Save predicted 2-page TIFF angle maps for each model.
% Each TIFF has:
%   page 1 = predicted fiber orientation angle
%   page 2 = predicted peptide-pitch angle
writePredictionTifs = true;

predRoot = fullfile(outRoot, 'Prediction_TIFs');
predDirPretrainedOnly = fullfile(predRoot, 'Pretrained_only_model');
predDirNoPretrain     = fullfile(predRoot, 'Without_pretrain_model');
predDirProposed       = fullfile(predRoot, 'Proposed_physics_pretrained_model');

if writePredictionTifs
    if ~isfolder(predRoot), mkdir(predRoot); end
    if ~isfolder(predDirPretrainedOnly), mkdir(predDirPretrainedOnly); end
    if ~isfolder(predDirNoPretrain), mkdir(predDirNoPretrain); end
    if ~isfolder(predDirProposed), mkdir(predDirProposed); end
end

% Exclude known unstable pixel-wise fitting case.
% Set {} to include all samples.
excludeSampleKeywords = {''};

threshold = 150;
fiberRange = 1:0.1:180;
pitchRange = 40:0.1:60;

% Reported peak values are snapped to the 0.1-degree histogram bin centers
% so values are consistent with the evaluation binning, e.g. 75.852 -> 75.85.
peakReportBinStep = 0.1;

% Optional debugging output.
writePerSampleCsv = true;
writeNumericSummaryCsv = true;

% Plot one representative ROI directly from the same prediction maps used
% for evaluation. This avoids mismatch between CSV and plot peaks.
plotRepresentativeFigure = true;

% Leave empty to plot the first ROI in fileList.
% Example: selectedPlotKeyword = 'Tendon_1_ROI01_TL';
selectedPlotKeyword = '';

% Save the generated figures as PNG and FIG in the output folder.
saveRepresentativeFigure = true;

%% ===================== MODEL PATHS: 64 experimental inputs =====================
pretrainedOnlyPath = 'model\SHG_Virtual_MNIST_Digit-IMAGE_Generic_Angle-ResNet-Filter-128-Epoch-200-Time-284624.mat';

noPrePath = 'model\TrainWithoutPretrain_64RealPatch-Epoch-200-Selected32Patch-64-TrainPatchSize-16-Balance-ColI50-ColII50-Seed-7-Time-260.mat';

proposedPath = 'model\Finetune_Proposed_64RealPatch-Epoch-200-Selected32Patch-64-TrainPatchSize-16-Balance-ColI50-ColII50-Seed-7-Time-260.mat';

%% ===================== CHECK INPUTS =====================
if ~isfolder(pshgFolder),  error('P-SHG folder does not exist: %s', pshgFolder); end
if ~isfolder(angleFolder), error('Angle folder does not exist: %s', angleFolder); end

fileList = [dir(fullfile(pshgFolder, '*.tif')); dir(fullfile(pshgFolder, '*.tiff'))];
if isempty(fileList), error('No TIFF files found in: %s', pshgFolder); end

[~, sortIdx] = sort({fileList.name});
fileList = fileList(sortIdx);
fileList = filterFilesByKeyword(fileList, excludeSampleKeywords);

fprintf('P-SHG folder : %s\n', pshgFolder);
fprintf('Angle folder : %s\n', angleFolder);
fprintf('Output folder: %s\n', outRoot);
fprintf('Using %d ROI TIFF files after exclusion.\n', numel(fileList));
if ~isempty(excludeSampleKeywords)
    fprintf('Excluded keywords: %s\n', strjoin(excludeSampleKeywords, ', '));
end
fprintf('\n');

theoryTimeMap = readTheoreticalTimeLog(angleFolder);

%% ===================== LOAD MODELS =====================
fprintf('Loading pretrained-only model...\n');
load(pretrainedOnlyPath, 'net');
Net_Simulation = net;
clear net;
fprintf('Pretrained-only model loaded.\n\n');

fprintf('Loading non-pretrained transfer model, 64 experimental inputs...\n');
load(noPrePath, 'net');
Net_NoPretrain = net;
clear net;
fprintf('Non-pretrained model loaded.\n\n');

fprintf('Loading proposed physics-pretrained transfer model, 64 experimental inputs...\n');
load(proposedPath, 'net');
Net_Finetune = net;
clear net;
fprintf('Proposed model loaded.\n\n');

%% ===================== EVALUATION =====================
Result = makeEmptyResultTable();

repData = struct();
repData.isSet = false;

for m = 1:numel(fileList)
    fileName = fileList(m).name;
    [~, sampleName, ext] = fileparts(fileName);

    pshgPath  = fullfile(pshgFolder, fileName);
    anglePath = fullfile(angleFolder, fileName);

    if ~isfile(anglePath)
        altExt = '.tif';
        if strcmpi(ext, '.tif'), altExt = '.tiff'; end

        altAnglePath = fullfile(angleFolder, [sampleName, altExt]);
        if isfile(altAnglePath)
            anglePath = altAnglePath;
        else
            warning('Missing angle TIFF for %s. Skipped.', sampleName);
            continue;
        end
    end

    fprintf('[%02d/%02d] %s\n', m, numel(fileList), fileName);

    try
        SHG_Image = double(tiffreadVolume(pshgPath));
        if size(SHG_Image,3) < 18
            error('Input stack has %d pages; expected at least 18.', size(SHG_Image,3));
        end
        SHG_Image = SHG_Image(:,:,1:18);

        infoA = imfinfo(anglePath);
        if numel(infoA) < 2
            error('Angle TIFF has %d pages; expected 2 pages.', numel(infoA));
        end

        Fiber_Theory   = double(imread(anglePath, 'Index', 1));
        Peptide_Theory = double(imread(anglePath, 'Index', 2));

        mask = sum(SHG_Image, 3) <= threshold;

        Fiber_Theory(mask) = NaN;
        Peptide_Theory(mask) = NaN;
        Fiber_Theory(Fiber_Theory <= 0 | Fiber_Theory >= 180) = NaN;
        Peptide_Theory(Peptide_Theory <= 0) = NaN;

        [Fiber_Sim,   Peptide_Sim,   Time_Sim]   = predict_image(Net_Simulation, SHG_Image, threshold);
        [Fiber_NoPre, Peptide_NoPre, Time_NoPre] = predict_image(Net_NoPretrain, SHG_Image, threshold);
        [Fiber_Fine,  Peptide_Fine,  Time_Fine]  = predict_image(Net_Finetune, SHG_Image, threshold);

        [Fiber_Sim,   Peptide_Sim]   = cleanPrediction(Fiber_Sim,   Peptide_Sim,   mask);
        [Fiber_NoPre, Peptide_NoPre] = cleanPrediction(Fiber_NoPre, Peptide_NoPre, mask);
        [Fiber_Fine,  Peptide_Fine]  = cleanPrediction(Fiber_Fine,  Peptide_Fine,  mask);

        % Store one representative ROI for plotting.
        % This uses the exact same cleaned maps that are passed to appendResultRow,
        % so the histogram peak in the figure matches the CSV evaluation.
        if plotRepresentativeFigure && ~repData.isSet
            if isempty(selectedPlotKeyword) || ~isempty(strfind(lower(fileName), lower(selectedPlotKeyword))) %#ok<STREMP>
                SHG_Image_ForPlot = SHG_Image;
                SHG_Image_ForPlot(mask) = 0;

                repData.isSet = true;
                repData.fileName = fileName;
                repData.SHG_Image = SHG_Image_ForPlot;
                repData.Fiber_Theory = Fiber_Theory;
                repData.Fiber_Sim = Fiber_Sim;
                repData.Fiber_NoPre = Fiber_NoPre;
                repData.Fiber_Fine = Fiber_Fine;
                repData.Peptide_Theory = Peptide_Theory;
                repData.Peptide_Sim = Peptide_Sim;
                repData.Peptide_NoPre = Peptide_NoPre;
                repData.Peptide_Fine = Peptide_Fine;

                % Store exact peak values using the same functions/ranges as CSV.
                % The plot will use these values for red lines/titles, so it
                % cannot drift from tendon_eval_per_sample.csv.
                repData.FiberPeaks = [ ...
                    snapPeakToBinCenter(histogramPeakCircular180(Fiber_Theory, fiberRange), peakReportBinStep), ...
                    snapPeakToBinCenter(histogramPeakCircular180(Fiber_Sim, fiberRange), peakReportBinStep), ...
                    snapPeakToBinCenter(histogramPeakCircular180(Fiber_NoPre, fiberRange), peakReportBinStep), ...
                    snapPeakToBinCenter(histogramPeakCircular180(Fiber_Fine, fiberRange), peakReportBinStep)];

                repData.PitchPeaks = [ ...
                    snapPeakToBinCenter(histogramPeak(Peptide_Theory, pitchRange), peakReportBinStep), ...
                    snapPeakToBinCenter(histogramPeak(Peptide_Sim, pitchRange), peakReportBinStep), ...
                    snapPeakToBinCenter(histogramPeak(Peptide_NoPre, pitchRange), peakReportBinStep), ...
                    snapPeakToBinCenter(histogramPeak(Peptide_Fine, pitchRange), peakReportBinStep)];
            end
        end

        % Save predicted angle maps as 2-page 32-bit float TIFFs.
        % Page 1 = fiber angle, page 2 = peptide-pitch angle.
        if writePredictionTifs
            writeTwoPageAngleTiff(fullfile(predDirPretrainedOnly, fileName), Fiber_Sim,   Peptide_Sim);
            writeTwoPageAngleTiff(fullfile(predDirNoPretrain,     fileName), Fiber_NoPre, Peptide_NoPre);
            writeTwoPageAngleTiff(fullfile(predDirProposed,       fileName), Fiber_Fine,  Peptide_Fine);
        end

        theoryTime = getTheoryTimeForSample(theoryTimeMap, sampleName);

        Result = appendResultRow(Result, sampleName, ...
            'Pixel-wise theoretical fitting', 'Conventional fitting', theoryTime, ...
            Fiber_Theory, Peptide_Theory, fiberRange, pitchRange, peakReportBinStep);

        Result = appendResultRow(Result, sampleName, ...
            'Pretrained-only model', 'Synthetic pretraining only', Time_Sim, ...
            Fiber_Sim, Peptide_Sim, fiberRange, pitchRange, peakReportBinStep);

        Result = appendResultRow(Result, sampleName, ...
            'Non-pretrained transfer model', '64 experimental inputs only', Time_NoPre, ...
            Fiber_NoPre, Peptide_NoPre, fiberRange, pitchRange, peakReportBinStep);

        Result = appendResultRow(Result, sampleName, ...
            'Proposed physics-pretrained transfer model', 'Synthetic pretraining + 64 experimental inputs', Time_Fine, ...
            Fiber_Fine, Peptide_Fine, fiberRange, pitchRange, peakReportBinStep);

    catch
        warning('Failed %s. Please check this file manually.', fileName);
    end
end

%% ===================== SUMMARY + PAPER TABLE =====================
Summary = buildSummaryTable(Result);
PaperTable = buildPaperTable(Summary);

paperPath = fullfile(outRoot, 'tendon_eval_paper_table.csv');
writetable(PaperTable, paperPath);
fprintf('\nSaved paper table: %s\n', paperPath);

if writeNumericSummaryCsv
    summaryPath = fullfile(outRoot, 'tendon_eval_numeric_summary.csv');

    % Save numeric summary with fixed 2 decimals, including trailing zeros.
    SummaryOut = formatTableNumeric2Decimals(Summary, {'N'});
    writetable(SummaryOut, summaryPath);

    fprintf('Saved numeric summary: %s\n', summaryPath);
end

if writePerSampleCsv
    perPath = fullfile(outRoot, 'tendon_eval_per_sample.csv');

    % Save per-sample values with fixed 2 decimals, including trailing zeros.
    ResultOut = formatTableNumeric2Decimals(Result, {'SelectedExperimentalInputs'});

    try
        writetable(ResultOut, perPath);
        fprintf('Saved per-sample table: %s\n', perPath);
    catch
        warning('Could not overwrite tendon_eval_per_sample.csv. The file may be open or locked.');
        backupPath = fullfile(outRoot, ['tendon_eval_per_sample_' datestr(now,'yyyymmdd_HHMMSS') '.csv']);
        writetable(ResultOut, backupPath);
        fprintf('Saved per-sample table to backup file instead: %s\n', backupPath);
    end
end

fprintf('\nDone.\n');

%% ===================== REPRESENTATIVE PLOT =====================
if plotRepresentativeFigure
    if repData.isSet
        fprintf('Plotting representative ROI: %s\n', repData.fileName);

        low_trim = 0;
        high_trim = 1;
        pLow = 5;
        pHigh = 90;
        RGB_scale = 0.3;

        % Fiber overlay display range.
        low = 1;
        high = 90;
        scale = 0.5;

        % Same rotation used in the separate Tendon_results_plot script.
        Rotate_Angle = 95;

        % Peak label in histograms is controlled inside plot_eval_synced_histogram:
        %   peakTextFontSize, peakLineWidth, peakTextYOffset.


        plot_sub_tendon_results( ...
            repData.SHG_Image, ...
            repData.Fiber_Theory, ...
            repData.Fiber_Sim, ...
            repData.Fiber_NoPre, ...
            repData.Fiber_Fine, ...
            repData.Peptide_Theory, ...
            repData.Peptide_Sim, ...
            repData.Peptide_NoPre, ...
            repData.Peptide_Fine, ...
            low_trim, ...
            high_trim, ...
            pLow, ...
            pHigh, ...
            RGB_scale, ...
            low, ...
            high, ...
            scale, ...
            Rotate_Angle, ...
            repData.FiberPeaks, ...
            repData.PitchPeaks ...
        );

        if saveRepresentativeFigure
            figOutDir = fullfile(outRoot, 'Representative_Plots');
            if ~isfolder(figOutDir), mkdir(figOutDir); end

            safeName = regexprep(repData.fileName, '[^\w\-.]', '_');
            safeName = regexprep(safeName, '\.tiff?$', '');

            figHandles = findall(0, 'Type', 'figure');
            figHandles = flipud(figHandles(:));

            for ff = 1:numel(figHandles)
                pngPath = fullfile(figOutDir, sprintf('%s_Figure_%02d.png', safeName, ff));
                figPath = fullfile(figOutDir, sprintf('%s_Figure_%02d.fig', safeName, ff));

                try
                    saveas(figHandles(ff), pngPath);
                    saveas(figHandles(ff), figPath);
                catch
                    warning('Could not save representative figure %d.', ff);
                end
            end

            fprintf('Saved representative plots to: %s\n', figOutDir);
        end
    else
        warning('No representative ROI was stored for plotting.');
    end
end

%% ===================== LOCAL FUNCTIONS =====================






function writeTwoPageAngleTiff(outPath, fiberMap, pitchMap)
    % Write a 2-page 32-bit float TIFF:
    %   page 1 = fiber orientation angle
    %   page 2 = peptide-pitch angle
    %
    % This version uses function-style Tiff calls for older MATLAB versions:
    %   setTag(t, tagstruct), write(t, img), writeDirectory(t), close(t)

    fiberMap = single(gather(fiberMap));
    pitchMap = single(gather(pitchMap));

    if isfile(outPath)
        delete(outPath);
    end

    t = Tiff(outPath, 'w');

    % Page 1: fiber
    tagstruct = makeFloatTiffTag(fiberMap);
    setTag(t, tagstruct);
    write(t, fiberMap);

    % Page 2: pitch
    writeDirectory(t);
    tagstruct = makeFloatTiffTag(pitchMap);
    setTag(t, tagstruct);
    write(t, pitchMap);

    close(t);

    % Quick check
    info = imfinfo(outPath);
    if numel(info) ~= 2
        warning('Expected 2 pages but wrote %d page(s): %s', numel(info), outPath);
    end
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


function T = formatTableNumeric2Decimals(T, keepNumericVars)
    % Convert numeric output columns to fixed 2-decimal strings.
    % This preserves trailing zeros in CSV, e.g. 70 -> 70.00.
    %
    % keepNumericVars are kept numeric/integer, e.g. N or SelectedExperimentalInputs.

    if nargin < 2
        keepNumericVars = {};
    end

    varNames = T.Properties.VariableNames;

    for i = 1:numel(varNames)
        v = varNames{i};

        if any(strcmp(v, keepNumericVars))
            continue;
        end

        if isnumeric(T.(v)) || islogical(T.(v))
            vals = double(gather(T.(v)));
            out = strings(size(vals));

            for j = 1:numel(vals)
                if isfinite(vals(j))
                    out(j) = string(sprintf('%.2f', vals(j)));
                else
                    out(j) = "NaN";
                end
            end

            T.(v) = out;
        end
    end
end


function fileList = filterFilesByKeyword(fileList, excludeKeywords)
    if isempty(excludeKeywords), return; end

    keep = true(numel(fileList), 1);
    for i = 1:numel(fileList)
        nm = lower(fileList(i).name);
        for k = 1:numel(excludeKeywords)
            key = lower(excludeKeywords{k});
            if ~isempty(strfind(nm, key)) %#ok<STREMP>
                keep(i) = false;
            end
        end
    end
    fileList = fileList(keep);
end

function T = makeEmptyResultTable()
    T = table('Size', [0 9], ...
        'VariableTypes', {'string','string','string','double','double','double','double','double','double'}, ...
        'VariableNames', {'Sample','Method','TrainingCondition','Time_sec', ...
                          'FiberPeak_deg','FiberTotalVariation', ...
                          'PeptidePitchPeak_deg','PitchTotalVariation', ...
                          'SelectedExperimentalInputs'});
end

function [fiberMap, pitchMap] = cleanPrediction(fiberMap, pitchMap, mask)
    fiberMap(mask) = NaN;
    pitchMap(mask) = NaN;

    fiberMap(fiberMap <= 0 | fiberMap >= 180) = NaN;
    pitchMap(pitchMap <= 0) = NaN;
end

function T = appendResultRow(T, sample, method, trainingCondition, timeSec, fiberMap, pitchMap, fiberRange, pitchRange, peakReportBinStep)
    fiberPeak = histogramPeakCircular180(fiberMap, fiberRange);
    pitchPeak = histogramPeak(pitchMap, pitchRange);

    % Snap reported peaks to histogram bin centers for consistent CSV/plot values.
    fiberPeak = snapPeakToBinCenter(fiberPeak, peakReportBinStep);
    pitchPeak = snapPeakToBinCenter(pitchPeak, peakReportBinStep);

    normFiber = normalizeImageForMetric(fiberMap);
    normPitch = normalizeImageForMetric(pitchMap);

    fiberTV = double(gather(computeTotalVariation(normFiber)));
    pitchTV = double(gather(computeTotalVariation(normPitch)));

    if strcmp(method, 'Pixel-wise theoretical fitting') || strcmp(method, 'Pretrained-only model')
        selectedInputs = 0;
    else
        selectedInputs = 64;
    end

    T = [T; {string(sample), string(method), string(trainingCondition), double(timeSec), ...
        fiberPeak, fiberTV, pitchPeak, pitchTV, selectedInputs}]; %#ok<AGROW>
end


function peakVal = snapPeakToBinCenter(peakVal, binStep)
    % Snap a continuous/fallback peak to the nearest histogram bin center.
    % For binStep = 0.1, valid reported centers are ... 75.750, 75.850, etc.
    if ~isfinite(peakVal)
        return;
    end

    peakVal = round((peakVal - binStep/2) ./ binStep) .* binStep + binStep/2;
    peakVal = double(peakVal);
end


function Summary = buildSummaryTable(Result)
    methodOrder = ["Pixel-wise theoretical fitting"; ...
                   "Pretrained-only model"; ...
                   "Non-pretrained transfer model"; ...
                   "Proposed physics-pretrained transfer model"];

    trainingMap = containers.Map( ...
        {'Pixel-wise theoretical fitting', 'Pretrained-only model', 'Non-pretrained transfer model', 'Proposed physics-pretrained transfer model'}, ...
        {'Conventional fitting', 'Synthetic pretraining only', '64 experimental inputs only', 'Synthetic pretraining + 64 experimental inputs'});

    Summary = table();

    for i = 1:numel(methodOrder)
        method = methodOrder(i);
        idx = Result.Method == method;

        if ~any(idx)
            continue;
        end

        row = table();
        row.Method = method;
        row.TrainingCondition = string(trainingMap(char(method)));
        row.N = sum(idx);

        row.Time_sec_Mean = mean(Result.Time_sec(idx), 'omitnan');
        row.Time_sec_SD   = std(Result.Time_sec(idx),  'omitnan');

        row.FiberPeak_deg_Mean = mean(Result.FiberPeak_deg(idx), 'omitnan');
        row.FiberPeak_deg_SD   = std(Result.FiberPeak_deg(idx),  'omitnan');

        row.FiberTotalVariation_Mean = mean(Result.FiberTotalVariation(idx), 'omitnan');
        row.FiberTotalVariation_SD   = std(Result.FiberTotalVariation(idx),  'omitnan');

        row.PeptidePitchPeak_deg_Mean = mean(Result.PeptidePitchPeak_deg(idx), 'omitnan');
        row.PeptidePitchPeak_deg_SD   = std(Result.PeptidePitchPeak_deg(idx),  'omitnan');

        row.PitchTotalVariation_Mean = mean(Result.PitchTotalVariation(idx), 'omitnan');
        row.PitchTotalVariation_SD   = std(Result.PitchTotalVariation(idx),  'omitnan');

        Summary = [Summary; row]; %#ok<AGROW>
    end

    % TV ratio = method TV mean / pixel-wise theoretical fitting TV mean.
    refIdx = Summary.Method == "Pixel-wise theoretical fitting";
    refFiberTV = Summary.FiberTotalVariation_Mean(refIdx);
    refPitchTV = Summary.PitchTotalVariation_Mean(refIdx);

    Summary.FiberTotalVariation_Mean = double(gather(Summary.FiberTotalVariation_Mean));
    Summary.FiberTotalVariation_SD   = double(gather(Summary.FiberTotalVariation_SD));
    Summary.PitchTotalVariation_Mean = double(gather(Summary.PitchTotalVariation_Mean));
    Summary.PitchTotalVariation_SD   = double(gather(Summary.PitchTotalVariation_SD));

    refFiberTV = double(gather(refFiberTV));
    refPitchTV = double(gather(refPitchTV));

    Summary.FiberTVRatio = double(gather(Summary.FiberTotalVariation_Mean ./ refFiberTV));
    Summary.PitchTVRatio = double(gather(Summary.PitchTotalVariation_Mean ./ refPitchTV));
end

function PaperTable = buildPaperTable(Summary)
    PaperTable = table();
    PaperTable.Method = Summary.Method;
    PaperTable.("Training condition") = Summary.TrainingCondition;

    PaperTable.("Time per ROI (s)") = composeMeanSD(Summary.Time_sec_Mean, Summary.Time_sec_SD, 2);
    PaperTable.("Fiber orientation peak (°)") = composeMeanSD(Summary.FiberPeak_deg_Mean, Summary.FiberPeak_deg_SD, 2);
    PaperTable.("Fiber orientation TV") = composeMeanSD(Summary.FiberTotalVariation_Mean, Summary.FiberTotalVariation_SD, 2);
    PaperTable.("Fiber orientation TV ratio") = compose('%.2f', double(gather(Summary.FiberTVRatio)));

    PaperTable.("Peptide-pitch angle peak (°)") = composeMeanSD(Summary.PeptidePitchPeak_deg_Mean, Summary.PeptidePitchPeak_deg_SD, 2);
    PaperTable.("Peptide-pitch angle TV") = composeMeanSD(Summary.PitchTotalVariation_Mean, Summary.PitchTotalVariation_SD, 2);
    PaperTable.("Peptide-pitch angle TV ratio") = compose('%.2f', double(gather(Summary.PitchTVRatio)));
end

function out = composeMeanSD(mu, sd, ndigits)
    mu = double(gather(mu));
    sd = double(gather(sd));
    fmt = sprintf('%%.%df ± %%.%df', ndigits, ndigits);
    out = strings(numel(mu), 1);
    for i = 1:numel(mu)
        if isfinite(mu(i)) && isfinite(sd(i))
            out(i) = string(sprintf(fmt, mu(i), sd(i)));
        elseif isfinite(mu(i))
            out(i) = string(sprintf(sprintf('%%.%df', ndigits), mu(i)));
        else
            out(i) = "NaN";
        end
    end
end

function peakVal = histogramPeakCircular180(img, rangeVals)
    %#ok<INUSD>
    x = double(img(:));
    x = x(isfinite(x));
    x = mod(x, 180);
    x = x(x > 0 & x < 180);

    if isempty(x)
        peakVal = NaN;
        return;
    end

    try
        edges0 = 0:0.1:180;
        counts0 = double(histcounts(x, edges0));
        centers0 = (edges0(1:end-1) + edges0(2:end)) ./ 2;
        [~, idx0] = max(counts0);
        roughPeak = centers0(idx0);

        shiftVal = 90 - roughPeak;
        xShift = mod(x + shiftVal, 180);
        xShift = xShift(xShift > 0 & xShift < 180);

        [~, ~, centerShift, ~] = fitGaussianHistogram(xShift, 0:1:180, 'angle', 0, []);
        if ~isfinite(centerShift)
            error('Circular Gaussian center is not finite.');
        end

        peakVal = mod(centerShift - shiftVal, 180);
        if peakVal < 0.5
            peakVal = 180;
        end
    catch
        edges0 = 0:0.1:180;
        counts0 = double(histcounts(x, edges0));
        centers0 = (edges0(1:end-1) + edges0(2:end)) ./ 2;
        [~, idx0] = max(counts0);
        peakVal = centers0(idx0);
    end
end

function peakVal = histogramPeak(img, rangeVals)
    x = double(img(:));
    x = x(isfinite(x));
    x = x(x >= min(rangeVals) & x <= max(rangeVals));

    if isempty(x)
        peakVal = NaN;
        return;
    end

    try
        [~, ~, peakVal, ~] = fitGaussianHistogram(x, rangeVals, 'angle', 0, []);
        if ~isfinite(peakVal)
            error('Gaussian center is not finite.');
        end
    catch
        binStep = median(diff(rangeVals));
        edges = [rangeVals(:); max(rangeVals) + binStep];
        counts = histcounts(x, edges);
        centers = (edges(1:end-1) + edges(2:end)) ./ 2;
        [~, idx] = max(counts);
        peakVal = centers(idx);
    end
end

function out = normalizeImageForMetric(img)
    img = double(img);
    valid = isfinite(img);

    if ~any(valid(:))
        out = zeros(size(img));
        return;
    end

    minVal = min(img(valid));
    maxVal = max(img(valid));

    out = (img - minVal) ./ (maxVal - minVal + eps);
    out(~valid) = 0;
end

function val = computeTotalVariation(img)
    dx = diff(img, 1, 2);
    dy = diff(img, 1, 1);

    dx = dx(1:end-1, :);
    dy = dy(:, 1:end-1);

    val = sum(sqrt(dx.^2 + dy.^2), 'all');
end

function timeMap = readTheoreticalTimeLog(angleFolder)
    timeMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    logPath = fullfile(angleFolder, 'theoretical_model_time_log.csv');

    if ~isfile(logPath)
        return;
    end

    try
        L = readtable(logPath, 'TextType', 'string');
        if all(ismember({'Sample','Time_sec'}, L.Properties.VariableNames))
            for i = 1:height(L)
                timeMap(char(L.Sample(i))) = double(L.Time_sec(i));
            end
        end
    catch
    end
end

function t = getTheoryTimeForSample(timeMap, sampleName)
    t = NaN;
    key = char(sampleName);

    if isKey(timeMap, key)
        t = timeMap(key);
    end
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

%% ===================== LOCAL PLOTTING FUNCTIONS =====================
function plot_sub_tendon_results(I_SHG, ...
                      Fiber_Angle_GT, ...
                      Fiber_Angle_Predict_Simulation, ...
                      Fiber_Angle_Predict_NoPretrain, ...
                      Fiber_Angle_Predict_Finetune, ...
                      Peptide_Angle_GT, ...
                      Peptide_Angle_Predict_Simulation, ...
                      Peptide_Angle_Predict_NoPretrain, ...
                      Peptide_Angle_Predict_Finetune, ...
                      low_trim, high_trim, pLow, pHigh, RGB_scale, ... %#ok<INUSD>
                      low, high, scale, Rotate_Angle, evalFiberPeaks, evalPitchPeaks)
% plot_sub_tendon_results_split_fiber_pitch
% Create TWO separate figures:
%   Figure 1: Fiber results  (4 rows x 3 columns)
%   Figure 2: Pitch results  (4 rows x 3 columns)
%
% Rows:
%   1 = Theoretical fitting
%   2 = Pretrained-only
%   3 = Without pretrain
%   4 = Proposed
%
% Columns for each figure:
%   1 = map
%   2 = surface
%   3 = peak / histogram

    if nargin < 19
        evalFiberPeaks = [];
    end
    if nargin < 20
        evalPitchPeaks = [];
    end

    %% ===================== Convert / rotate =====================
    Fiber_Angle_GT                 = double(gather(Fiber_Angle_GT));
    Fiber_Angle_Predict_Simulation = double(gather(Fiber_Angle_Predict_Simulation));
    Fiber_Angle_Predict_NoPretrain = double(gather(Fiber_Angle_Predict_NoPretrain));
    Fiber_Angle_Predict_Finetune   = double(gather(Fiber_Angle_Predict_Finetune));

    Peptide_Angle_GT                 = double(gather(Peptide_Angle_GT));
    Peptide_Angle_Predict_Simulation = double(gather(Peptide_Angle_Predict_Simulation));
    Peptide_Angle_Predict_NoPretrain = double(gather(Peptide_Angle_Predict_NoPretrain));
    Peptide_Angle_Predict_Finetune   = double(gather(Peptide_Angle_Predict_Finetune));

    I_SHG_Rot = rotate_and_crop_image(I_SHG, Rotate_Angle);

    %% ===================== Match run_tendon_eval mask for histogram =====================
    % The show script passes SHG_Image after background masking, where masked pixels are 0.
    % For histogram/peak evaluation, remove those same background pixels from all maps.
    histMask = compute_hist_mask_from_pshg(I_SHG);

    Fiber_Angle_GT(histMask)                 = NaN;
    Fiber_Angle_Predict_Simulation(histMask) = NaN;
    Fiber_Angle_Predict_NoPretrain(histMask) = NaN;
    Fiber_Angle_Predict_Finetune(histMask)   = NaN;

    Peptide_Angle_GT(histMask)                 = NaN;
    Peptide_Angle_Predict_Simulation(histMask) = NaN;
    Peptide_Angle_Predict_NoPretrain(histMask) = NaN;
    Peptide_Angle_Predict_Finetune(histMask)   = NaN;

    Fiber_Angle_GT(Fiber_Angle_GT <= 0 | Fiber_Angle_GT >= 180) = NaN;
    Fiber_Angle_Predict_Simulation(Fiber_Angle_Predict_Simulation <= 0 | Fiber_Angle_Predict_Simulation >= 180) = NaN;
    Fiber_Angle_Predict_NoPretrain(Fiber_Angle_Predict_NoPretrain <= 0 | Fiber_Angle_Predict_NoPretrain >= 180) = NaN;
    Fiber_Angle_Predict_Finetune(Fiber_Angle_Predict_Finetune <= 0 | Fiber_Angle_Predict_Finetune >= 180) = NaN;

    Peptide_Angle_GT(Peptide_Angle_GT <= 0) = NaN;
    Peptide_Angle_Predict_Simulation(Peptide_Angle_Predict_Simulation <= 0) = NaN;
    Peptide_Angle_Predict_NoPretrain(Peptide_Angle_Predict_NoPretrain <= 0) = NaN;
    Peptide_Angle_Predict_Finetune(Peptide_Angle_Predict_Finetune <= 0) = NaN;


    Fiber_GT_Rot    = rotate_and_crop_image(Fiber_Angle_GT, Rotate_Angle);
    Fiber_Sim_Rot   = rotate_and_crop_image(Fiber_Angle_Predict_Simulation, Rotate_Angle);
    Fiber_NoPre_Rot = rotate_and_crop_image(Fiber_Angle_Predict_NoPretrain, Rotate_Angle);
    Fiber_Fine_Rot  = rotate_and_crop_image(Fiber_Angle_Predict_Finetune, Rotate_Angle);

    Pitch_GT_Rot    = rotate_and_crop_image(Peptide_Angle_GT, Rotate_Angle);
    Pitch_Sim_Rot   = rotate_and_crop_image(Peptide_Angle_Predict_Simulation, Rotate_Angle);
    Pitch_NoPre_Rot = rotate_and_crop_image(Peptide_Angle_Predict_NoPretrain, Rotate_Angle);
    Pitch_Fine_Rot  = rotate_and_crop_image(Peptide_Angle_Predict_Finetune, Rotate_Angle);

    %% ===================== Groups =====================
    fiber_maps = {Fiber_GT_Rot, Fiber_Sim_Rot, Fiber_NoPre_Rot, Fiber_Fine_Rot};
    fiber_hist = {Fiber_Angle_GT, Fiber_Angle_Predict_Simulation, Fiber_Angle_Predict_NoPretrain, Fiber_Angle_Predict_Finetune};

    pitch_maps = {Pitch_GT_Rot, Pitch_Sim_Rot, Pitch_NoPre_Rot, Pitch_Fine_Rot};
    pitch_hist = {Peptide_Angle_GT, Peptide_Angle_Predict_Simulation, Peptide_Angle_Predict_NoPretrain, Peptide_Angle_Predict_Finetune};

    row_names = {'Theoretical fitting', 'Pretrained-only', 'Without pretrain', 'Proposed'};

    cmap_sym   = custom_colormap(256);
    map_fiber  = cmap_sym;
    map_pitch  = cmap_sym;

    fiber_clim = [0 180];
    pitch_clim = [40 60];

    %% ===================== Fiber RGB maps =====================
    fiber_rgb = cell(1,4);
    for k = 1:4
        tmp = coat_SHG_with_angle(I_SHG_Rot, fiber_maps{k}, map_fiber, low, high, scale);
        fiber_rgb{k} = im2double(tmp(:,:,:,min(9, size(tmp,4))));
    end

    %% ===================== Histogram y-limits =====================
    yMax_fiber = get_common_hist_ymax(fiber_hist, 1:180);
    yMax_pitch = get_common_hist_ymax(pitch_hist, 40:60);

    %% ===================== FIGURE 0: P-SHG intensity =====================
    fig0 = figure('Units', 'pixels', ...
                  'Position', [20, 80, 520, 420], ...
                  'Color', 'w', ...
                  'Name', 'P-SHG intensity');

    intensityImg = flipud(compute_pshg_intensity(I_SHG_Rot));
    ax0 = axes('Parent', fig0, 'Position', [0.08, 0.14, 0.66, 0.76]);
    imagesc(ax0, intensityImg, [0 1]);
    axis(ax0, 'image', 'off');
    colormap(ax0, gray(256));
    set(ax0, 'Color', 'k');
    hold(ax0, 'on');
    rectangle(ax0, 'Position', [0.5, 0.5, size(intensityImg,2)-1, size(intensityImg,1)-1], ...
              'EdgeColor', 'w', 'LineWidth', 1.5, 'Clipping', 'on');
    hold(ax0, 'off');

    cb0 = colorbar(ax0, 'eastoutside');
    cb0.Position = [0.77, 0.14, 0.05, 0.76];
    cb0.Ticks = 0:0.2:1;
    cb0.TickLabels = {'0','0.2','0.4','0.6','0.8','1'};
    cb0.FontSize = 12;
    cb0.Label.String = 'P-SHG intensity';
    cb0.Label.FontSize = 18;
    cb0.Label.Rotation = 90;

    %% ===================== FIGURE 1: FIBER =====================
    fig1 = figure('Units', 'pixels', ...
                  'Position', [30, 40, 1480, 1120], ...
                  'Color', 'w', ...
                  'Name', 'Fiber comparison');

    row_h = 0.158;
    row_gap = 0.060;
    y0 = 0.070;
    y_pos = y0 + (3:-1:0) * (row_h + row_gap);

    x1 = 0.055; w1 = 0.165;  % map
    x2 = 0.275; w2 = 0.250;  % surface
    cbx = 0.540; cbw = 0.012;
    x3 = 0.595; w3 = 0.205;  % histogram
    cby = y0;
    cbh = 4*row_h + 3*row_gap;

    for r = 1:4
        y = y_pos(r);

        ax1 = axes('Parent', fig1, 'Position', [x1, y, w1, row_h]);
        plot_fiber_rgb_no_roi(ax1, fiber_rgb{r});

        ax2 = axes('Parent', fig1, 'Position', [x2, y, w2, row_h]);
        plot_angle_surf_full_rgb(ax2, fiber_maps{r}, fiber_rgb{r}, fiber_clim, '\theta (deg)');

        ax3 = axes('Parent', fig1, 'Position', [x3, y, w3, row_h]);
        plot_eval_synced_histogram(ax3, fiber_hist{r}, 1:180, 'Fiber orientation (\circ)', yMax_fiber, 'fiber', get_peak_from_vector(evalFiberPeaks, r));
        set(ax3, 'FontSize', 9, 'LineWidth', 1);
    end

    cb_ax1 = axes('Parent', fig1, 'Position', [cbx, cby, cbw, cbh], 'Visible', 'off');
    colormap(cb_ax1, map_fiber); caxis(cb_ax1, fiber_clim);
    cb1 = colorbar(cb_ax1, 'Location', 'eastoutside');
    cb1.Units = 'normalized';
    cb1.Position = [cbx, cby, cbw, cbh];
    cb1.Ticks = 0:30:180;
    cb1.FontSize = 10;
    cb1.Label.String = 'Fiber orientation (\circ)';
    cb1.Label.FontSize = 11;

    %% ===================== FIGURE 2: PITCH =====================
    fig2 = figure('Units', 'pixels', ...
                  'Position', [80, 40, 1480, 1120], ...
                  'Color', 'w', ...
                  'Name', 'Pitch comparison');

    for r = 1:4
        y = y_pos(r);

        ax4 = axes('Parent', fig2, 'Position', [x1, y, w1, row_h]);
        plot_angle_map(ax4, pitch_maps{r}, map_pitch, pitch_clim);

        ax5 = axes('Parent', fig2, 'Position', [x2, y, w2, row_h]);
        plot_angle_surf_full_cmap(ax5, pitch_maps{r}, pitch_clim, '\alpha (deg)', map_pitch);

        ax6 = axes('Parent', fig2, 'Position', [x3, y, w3, row_h]);
        plot_eval_synced_histogram(ax6, pitch_hist{r}, 40:60, 'Peptide-pitch angle (\circ)', yMax_pitch, 'pitch', get_peak_from_vector(evalPitchPeaks, r));
        set(ax6, 'FontSize', 9, 'LineWidth', 1);
    end

    cb_ax2 = axes('Parent', fig2, 'Position', [cbx, cby, cbw, cbh], 'Visible', 'off');
    colormap(cb_ax2, map_pitch); caxis(cb_ax2, pitch_clim);
    cb2 = colorbar(cb_ax2, 'Location', 'eastoutside');
    cb2.Units = 'normalized';
    cb2.Position = [cbx, cby, cbw, cbh];
    cb2.Ticks = 40:5:60;
    cb2.FontSize = 10;
    cb2.Label.String = 'Peptide-pitch angle (\circ)';
    cb2.Label.FontSize = 11;
end

%% ========================= Helper functions =========================


function mask = compute_hist_mask_from_pshg(I_SHG)
    I_SHG = double(gather(I_SHG));

    if ndims(I_SHG) == 2
        mask = I_SHG <= 0;
    else
        mask = sum(I_SHG, 3) <= 0;
    end

    mask = logical(mask);
end

function val = get_peak_from_vector(v, idx)
    if isempty(v) || numel(v) < idx
        val = NaN;
    else
        val = v(idx);
    end
end

function plot_eval_synced_histogram(ax, img, rangeVals, xLabelStr, yMax, modeStr, peakOverride)
    % Histogram peak is computed using exactly the same functions as run_tendon_eval:
    %   Fiber: histogramPeakCircular180
    %   Pitch: histogramPeak
    %
    % Peak label uses axes-attached text placed above the histogram.
    % This avoids title/annotation copy-export problems.

    peakTextFontSize = 12;
    peakTextFontWeight = 'bold';
    peakLineWidth = 1.5;
    peakTextYOffset = 1.0;   % 1.04 closer, 1.12 farther above the plot

    axes(ax); %#ok<LAXES>
    x = clean_vec_for_hist(img, [min(rangeVals), max(rangeVals)]);

    if isempty(x)
        histogram(ax, NaN, rangeVals, 'EdgeColor', 'none', 'FaceAlpha', 1);
        peakVal = NaN;
    else
        histogram(ax, x, rangeVals, 'EdgeColor', 'none', 'FaceAlpha', 1);

        if nargin >= 7 && ~isempty(peakOverride) && isfinite(peakOverride)
            peakVal = peakOverride;
        else
            if strcmpi(modeStr, 'fiber')
                peakVal = histogramPeakCircular180(img, rangeVals);
            else
                peakVal = histogramPeak(img, rangeVals);
            end
        end
    end

    hold(ax, 'on');

    if isfinite(peakVal)
        xline(ax, peakVal, 'r', 'LineWidth', peakLineWidth);
        peakText = ['Peak = ' num2str(peakVal, '%.2f') '\circ'];
    else
        peakText = 'Peak = NaN';
    end

    % Do not use title() here; use axes-normalized text instead.
    title(ax, '');

    text(ax, 0.5, peakTextYOffset, peakText, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', peakTextFontSize, ...
        'FontWeight', peakTextFontWeight, ...
        'Color', 'k', ...
        'Interpreter', 'tex', ...
        'Clipping', 'off');

    xlabel(ax, xLabelStr);
    ylabel(ax, 'Counts');
    xlim(ax, [min(rangeVals), max(rangeVals)]);

    if nargin >= 5 && ~isempty(yMax) && isfinite(yMax) && yMax > 0
        ylim(ax, [0, yMax]);
    end

    axis(ax, 'square');
    hold(ax, 'off');
end


function yMax = get_common_hist_ymax(dataCells, edges)
    yMax = 1;
    for k = 1:numel(dataCells)
        x = clean_vec_for_hist(dataCells{k}, [min(edges), max(edges)]);
        if isempty(x), continue; end
        counts = histcounts(x, edges);
        yMax = max(yMax, max(counts));
    end
    yMax = yMax * 1.10;
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
              'EdgeColor', 'k', 'LineWidth', 1.0, 'Clipping', 'on');
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
              'EdgeColor', 'k', 'LineWidth', 1.0, 'Clipping', 'on');
    hold(ax, 'off');
end

function plot_angle_surf_full_rgb(ax, angleMap, rgbImg, climRange, zlabelStr)
    angleMap = double(gather(angleMap));
    angleMap(angleMap <= climRange(1) | angleMap >= climRange(2)) = NaN;

    [ny, nx] = size(angleMap);
    x_um = linspace(0, 200, nx);
    y_um = linspace(0, 200, ny);
    [X, Y] = meshgrid(x_um, y_um);

    rgbImg = im2double(rgbImg);
    h = surf(ax, X, Y, angleMap, rgbImg, 'EdgeColor', 'none');
    h.FaceColor = 'texturemap';

    zlim(ax, climRange);
    view(ax, 48, 58);
    set(ax, 'YDir', 'reverse');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'X (\mum)', 'FontSize', 8);
    ylabel(ax, 'Y (\mum)', 'FontSize', 8);
    zlabel(ax, zlabelStr, 'FontSize', 8);
    set(ax, 'FontSize', 8, 'LineWidth', 0.8);
end

function plot_angle_surf_full_cmap(ax, angleMap, climRange, zlabelStr, cmap)
    angleMap = double(gather(angleMap));
    angleMap(angleMap <= climRange(1) | angleMap >= climRange(2)) = NaN;

    [ny, nx] = size(angleMap);
    x_um = linspace(0, 200, nx);
    y_um = linspace(0, 200, ny);
    [X, Y] = meshgrid(x_um, y_um);

    surf(ax, X, Y, angleMap, angleMap, 'EdgeColor', 'none', 'FaceColor', 'interp');
    colormap(ax, cmap);
    caxis(ax, climRange);
    zlim(ax, climRange);
    view(ax, 48, 58);
    set(ax, 'YDir', 'reverse');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'X (\mum)', 'FontSize', 8);
    ylabel(ax, 'Y (\mum)', 'FontSize', 8);
    zlabel(ax, zlabelStr, 'FontSize', 8);
    set(ax, 'FontSize', 8, 'LineWidth', 0.8);
    shading(ax, 'interp');
end

function out = compute_pshg_intensity(Iin)
    Iin = double(gather(Iin));
    Iin(Iin < 0) = 0;

    if ndims(Iin) == 2
        out = Iin;
    else
        % For a P-SHG stack, use the mean intensity over all input states.
        out = mean(Iin, ndims(Iin));
        out = squeeze(out);
        while ndims(out) > 2
            out = mean(out, ndims(out));
        end
    end

    out = flipud(out);
    out = out - min(out(:));
    if max(out(:)) > 0
        out = out ./ max(out(:));
    end
end
