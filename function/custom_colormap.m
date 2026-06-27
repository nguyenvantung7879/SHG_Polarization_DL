function cmap = custom_colormap(n)
%PURPLE_RED_YELLOW_COLORMAP Generate a strong vivid colormap from purple to yellow.
%   cmap = purple_red_yellow_colormap(n) returns an n-by-3 matrix.

if nargin < 1
    n = 256; % default
end

% Define vivid RGB control points
controlColors = [...
    0.7 0.0 0.8 ;  % Strong Purple
    1.0 0.2 0.0 ;  % Pure Red
    1.0 0.4 0.0 ;  % Vivid Orange
    1.0 1.0 0.0;];  % Bright Yellow
controlColors = jet(256);
controlColors = flipud(controlColors);

% Control points normalized to [0,1]
controlPoints = linspace(0, 1, size(controlColors, 1));

% Interpolate across the full range
queryPoints = linspace(0, 1, n);

% Interpolate each RGB channel separately
cmap = interp1(controlPoints, controlColors, queryPoints, 'linear');
end
