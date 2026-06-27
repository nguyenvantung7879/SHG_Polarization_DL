function RGB_SHG = coat_SHG_with_angle(I_SHG, Angle_Map, colormapName, low, high, scale)
% COAT_SHG_WITH_ANGLE_CONTRAST_180
% Enhance SHG intensity contrast and coat it with a fiber angle color map.
%
% Inputs:
%   I_SHG        - SHG intensity volume, size [H x W x 18]
%   Angle_Map    - fiber angle map [H x W], values from 0 to 180 degrees
%   colormapName - name of colormap string (e.g., 'jet', 'hsv') or custom [N x 3]
%
% Output:
%   RGB_SHG      - Coated RGB volume [H x W x 3 x 18]
    [h, w, c] = size(I_SHG);  % Get size (H, W, 18)

    %% Step 1: Contrast enhance SHG using percentile clipping
    I_SHG = double(I_SHG);
    low_clip = prctile(I_SHG(:), low);
    high_clip = prctile(I_SHG(:), high);
    I_SHG_clipped = min(max(I_SHG, low_clip), high_clip);
    I_SHG_norm = (I_SHG_clipped - low_clip) / (high_clip - low_clip);

    %% Step 2: Normalize angle map from [0, 180]
    angle_norm = Angle_Map ./ 180;
    angle_norm = max(min(angle_norm, 1), 0);  % Clip to [0,1]

    % Convert angle to RGB using colormap
    if ischar(colormapName) || isstring(colormapName)
        cmap = colormap(colormapName);
    else
        cmap = colormapName;  % assume it's already [N x 3]
    end
    nColors = size(cmap, 1);

    % Map normalized angle to RGB
    idx = round(angle_norm * (nColors - 1)) + 1;
    RGB_angle = ind2rgb(idx, cmap);  % [H x W x 3]

    %% Step 3: Coat SHG slices with angle RGB
    RGB_SHG = zeros(h, w, 3, c);  % Output [H x W x 3 x 18]

    for k = 1:c
        for ch = 1:3
            RGB_SHG(:, :, ch, k) = I_SHG_norm(:, :, k) .* RGB_angle(:, :, ch);
        end
    end

    intensitySum = sum(I_SHG, 3);  % Sum over all channels
    backgroundMask = intensitySum < 150;  % Adjust threshold if needed
    
   % Loop through all slices to set background color
    for k = 1:c
        for ch = 1:3
            plane = RGB_SHG(:, :, ch, k);
            switch ch
                case 1, colorVal = 0;        % Red
                case 2, colorVal = 0;   % Green
                case 3, colorVal = 0;   % Blue
            end
            plane(backgroundMask) = colorVal;
            RGB_SHG(:, :, ch, k) = plane;
        end
    end
    RGB_SHG = RGB_SHG.^scale;
end
