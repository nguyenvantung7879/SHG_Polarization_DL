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
                      low, high, scale, Rotate_Angle)
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

    cmap_sym   = purple_red_yellow_colormap(256);
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

    intensityImg = compute_pshg_intensity(I_SHG_Rot);
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
        plot_eval_synced_histogram(ax3, fiber_hist{r}, 1:180, 'Fiber orientation (\circ)', yMax_fiber, 'fiber');
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
        plot_eval_synced_histogram(ax6, pitch_hist{r}, 40:60, 'Peptide-pitch angle (\circ)', yMax_pitch, 'pitch');
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

function plot_eval_synced_histogram(ax, img, rangeVals, xLabelStr, yMax, modeStr)
    % Histogram peak is computed using exactly the same functions as run_tendon_eval:
    %   Fiber: histogramPeakCircular180
    %   Pitch: histogramPeak

    axes(ax); %#ok<LAXES>
    x = clean_vec_for_hist(img, [min(rangeVals), max(rangeVals)]);

    if isempty(x)
        histogram(ax, NaN, rangeVals, 'EdgeColor', 'none', 'FaceAlpha', 1);
        peakVal = NaN;
    else
        histogram(ax, x, rangeVals, 'EdgeColor', 'none', 'FaceAlpha', 1);

        if strcmpi(modeStr, 'fiber')
            peakVal = histogramPeakCircular180(img, rangeVals);
        else
            peakVal = histogramPeak(img, rangeVals);
        end
    end

    hold(ax, 'on');

    if isfinite(peakVal)
        xline(ax, peakVal, 'r', 'LineWidth', 1.2);
        title(ax, ['Peak = ' num2str(peakVal, '%.3f') '\circ']);
    else
        title(ax, 'Peak = NaN');
    end

    xlabel(ax, xLabelStr);
    ylabel(ax, 'Counts');
    xlim(ax, [min(rangeVals), max(rangeVals)]);

    if nargin >= 5 && ~isempty(yMax) && isfinite(yMax) && yMax > 0
        ylim(ax, [0, yMax]);
    end

    axis(ax, 'square');
    hold(ax, 'off');
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
        edges0 = 0:1:180;
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
        edges0 = 0:1:180;
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
