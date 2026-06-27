%% generate_pitch_lookup_pq.m
% Generate lookup files required by generate_synthetic_pshg_mnist.m:
%   model/Pitch_Angle_Lookup.mat
%   model/P_Q_Value.mat
%
% Required helper:
%   Generic_SHG_Model.m
%
% This lookup maps generic P-SHG model parameters p and q to an equivalent
% peptide-pitch angle by fitting the generic response with the single-axis
% molecular model form.

clc;
clear;
close all;
format long e;
warning on;

%% ===================== PATH SETUP =====================
% Expected location:
%   repo/src/synthetic_data_generation/generate_pitch_lookup_pq.m
thisFile = mfilename('fullpath');
scriptFolder = fileparts(thisFile);
repoRoot = fileparts(fileparts(scriptFolder));   % repo root if script is in src/synthetic_data_generation

if isfolder(fullfile(repoRoot, 'src'))
    addpath(genpath(fullfile(repoRoot, 'src')));
end

if isfolder(fullfile(repoRoot, 'function'))
    addpath(genpath(fullfile(repoRoot, 'function')));
end

if exist('Generic_SHG_Model', 'file') ~= 2
    error('Missing required helper function: Generic_SHG_Model.m');
end

%% ===================== OUTPUT FOLDER =====================
modelFolder = fullfile(repoRoot, "SHG_Polarization_DL\data_simulation");
if ~isfolder(modelFolder)
    mkdir(modelFolder);
end

%% ===================== P/Q GRID =====================
% These ranges match the synthetic-data generation setting:
%   p in [-1, 0]
%   q in [-0.2, 0.2]
P = -(0:100) ./ 100;      % 0, -0.01, ..., -1
Q = (-20:20) ./ 100;      % -0.20, ..., 0.20

P = P(:);                 % column vector
Q = Q(:);                 % column vector

%% ===================== FIT MODEL =====================
% Fit generic P-SHG response with a molecular model-like expression.
% b is then converted to peptide-pitch angle using atan(2/b).
ft = fittype( ...
    'a*((sind(x).^2 + b .* cosd(x).^2).^2 + (2 .* sind(x) .* cosd(x)).^2)', ...
    'independent', 'x', ...
    'dependent', 'y');

opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.Lower = [0 0];
opts.StartPoint = [0.5 0.1];

x = 0:1:180;
Pitch_Angle_Table = NaN(numel(P), numel(Q));

fprintf('Generating Pitch_Angle_Lookup table...\n');
fprintf('P values: %d | Q values: %d | Total fits: %d\n', numel(P), numel(Q), numel(P)*numel(Q));

for m = 1:numel(P)
    if mod(m, 10) == 0 || m == 1 || m == numel(P)
        fprintf('P index %d/%d\n', m, numel(P));
    end

    for n = 1:numel(Q)
        I_SHG = Generic_SHG_Model(x, 0.0, P(m), Q(n));
        I_SHG = double(I_SHG(:)).';

        if max(I_SHG(:)) > 0
            I_SHG = I_SHG ./ max(I_SHG(:));
        else
            Pitch_Angle_Table(m, n) = NaN;
            continue;
        end

        try
            [xData, yData] = prepareCurveData(x, I_SHG);
            [fitresult, ~] = fit(xData, yData, ft, opts);
            PV = coeffvalues(fitresult);

            % Equivalent peptide-pitch angle.
            Pitch_Angle_Table(m, n) = atand(2 ./ PV(2));

        catch ME
            warning('Fit failed at P(%d)=%.3f, Q(%d)=%.3f: %s', ...
                m, P(m), n, Q(n), ME.message);
            Pitch_Angle_Table(m, n) = NaN;
        end
    end
end

%% ===================== SAVE LOOKUP FILES =====================
lookupPath = fullfile(modelFolder, 'Pitch_Angle_Lookup.mat');
pqPath     = fullfile(modelFolder, 'P_Q_Value.mat');

save(lookupPath, 'Pitch_Angle_Table', '-v7.3', '-nocompression');
save(pqPath, 'P', 'Q', '-v7.3', '-nocompression');

fprintf('\nDone.\n');
fprintf('Saved: %s\n', lookupPath);
fprintf('Saved: %s\n', pqPath);

%% ===================== QUICK CHECK =====================
fprintf('\nLookup table size: %d x %d\n', size(Pitch_Angle_Table,1), size(Pitch_Angle_Table,2));
fprintf('Pitch angle range: %.3f to %.3f deg\n', ...
    min(Pitch_Angle_Table(:), [], 'omitnan'), ...
    max(Pitch_Angle_Table(:), [], 'omitnan'));
