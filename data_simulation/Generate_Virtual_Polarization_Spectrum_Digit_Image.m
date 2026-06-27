%% generate_synthetic_pshg_mnist_generic_only.m
% Generate synthetic P-SHG image stacks using MNIST masks and the
% physics-based generic P-SHG model.
%
% Outputs:
%   SHG_Polarization_DL/data_simulation/Synthetic_PSHG_MNIST/
%       Polarization_Map/   SHG_IM_*.mat
%       Angle_Map/          Angle_IM_*.mat
%
% Each SHG_IM is H x W x 18.
% Each Angle_IM is H x W x 2:
%   channel 1 = collagen fiber orientation angle
%   channel 2 = peptide-pitch angle
%
% Required helper function:
%   Generic_SHG_Model.m
%
% Required lookup files:
%   model/Pitch_Angle_Lookup.mat
%   model/P_Q_Value.mat
%
% Note:
%   This generic-only version generates exactly N synthetic P-SHG stacks.
%   For N = 1.4e5, the output contains 140,000 polarization stacks and
%   140,000 corresponding angle-label maps.

clc;
clear;
close all;
format long e;
warning on;

%% ===================== PATH SETUP =====================
% This makes the script portable for GitHub.
% Expected location:
%   repo/src/synthetic_data_generation/generate_synthetic_pshg_mnist_generic_only.m

thisFile = mfilename('fullpath');
scriptFolder = fileparts(thisFile);

% If the file is stored in repo/src/synthetic_data_generation, this returns repo.
repoRoot = fileparts(fileparts(scriptFolder));

% If this script is run from another location, repoRoot may need to be changed manually.
if ~isfolder(fullfile(repoRoot, 'src')) && isfolder('src')
    repoRoot = pwd;
end

if isfolder(fullfile(repoRoot, 'src'))
    addpath(genpath(fullfile(repoRoot, 'src')));
end

% Optional function folder for older local structure.
if isfolder(fullfile(repoRoot, 'function'))
    addpath(genpath(fullfile(repoRoot, 'function')));
end

%% ===================== REPRODUCIBILITY =====================
% Fixed seed for reproducible GitHub results.
rng(7, 'twister');

%% ===================== LOAD MNIST MASK DATA =====================
digitDatasetPath = fullfile(matlabroot, 'toolbox', 'nnet', 'nndemos', ...
    'nndatasets', 'DigitDataset');

if ~isfolder(digitDatasetPath)
    error('MATLAB digit dataset not found: %s', digitDatasetPath);
end

imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders', true, ...
    'LabelSource', 'foldernames');

%% ===================== USER SETTINGS =====================
Save_fold = fullfile(repoRoot, 'SHG_Polarization_DL', 'data_simulation', 'Synthetic_PSHG_MNIST');

theta = 0:10:170;        % 18 excitation polarization angles, degrees
Dark_Counts = 2;         % dark noise counts
I_max = (5:5:4000).';    % intensity lookup as column vector

% Full dataset used in the manuscript.
N = 1.4e5;

% For quick GitHub/code testing, uncomment this:
% N = 200;
N = 1;

%% ===================== OUTPUT FOLDERS =====================
polarizationFolder = 'Polarization_Map';
angleFolder        = 'Angle_Map';

makeFolderIfMissing(Save_fold);
makeFolderIfMissing(fullfile(Save_fold, polarizationFolder));
makeFolderIfMissing(fullfile(Save_fold, angleFolder));

%% ===================== CHECK REQUIRED FILES/FUNCTIONS =====================
if exist('Generic_SHG_Model', 'file') ~= 2
    error('Missing required helper function: Generic_SHG_Model.m');
end

lookupPath = fullfile(repoRoot, 'model', 'Pitch_Angle_Lookup.mat');
pqPath     = fullfile(repoRoot, 'model', 'P_Q_Value.mat');

% Fallback: allow lookup files in current working directory.
if ~isfile(lookupPath) && isfile('Pitch_Angle_Lookup.mat')
    lookupPath = 'Pitch_Angle_Lookup.mat';
end
if ~isfile(pqPath) && isfile('P_Q_Value.mat')
    pqPath = 'P_Q_Value.mat';
end

if ~isfile(lookupPath)
    error('Missing lookup file: %s', lookupPath);
end
if ~isfile(pqPath)
    error('Missing P/Q file: %s', pqPath);
end

load(lookupPath, 'Pitch_Angle_Table');
load(pqPath, 'P', 'Q');

P = P(:);
Q = Q(:);

%% ===================== READ MNIST MASK TEMPLATES =====================
file_num = numel(imds.Files);
if file_num == 0
    error('No MNIST files found in: %s', digitDatasetPath);
end

fprintf('Reading %d MNIST mask templates...\n', file_num);

for m = 1:file_num
    file_name = imds.Files{m};
    Temp = imread(file_name);

    if m == 1
        NX = size(Temp, 1);
        NY = size(Temp, 2);
        IM = zeros(NX, NY, file_num, 'double');
    end

    Mask = zeros(NX, NY, 'double');
    Mask(Temp > 10) = 1;
    IM(:, :, m) = Mask;
end

%% ===================== BUILD RANDOM TEMPLATE INDEX LIST =====================
NP = fix(N / file_num);
IND = [];

for m = 1:NP
    idx = randperm(file_num);
    IND = [IND, idx]; %#ok<AGROW>
end

Residual = mod(N, file_num);
if Residual > 0
    idx = randperm(file_num);
    IND = [IND, idx(1:Residual)]; %#ok<AGROW>
end

IND = IND(:);
clear NP idx Residual;

%% ===================== GENERATE GENERIC SYNTHETIC P-SHG DATA =====================
N_theta = numel(theta);
NI = numel(I_max);

fprintf('Generating %d generic synthetic P-SHG stacks...\n', N);
fprintf('Output folder: %s\n\n', Save_fold);

for m = 1:N
    if mod(m, 1000) == 0 || m == 1 || m == N
        fprintf('[%d/%d]\n', m, N);
    end

    Mask = IM(:, :, IND(m));
    NL = sum(Mask, "all");

    if NL == 0
        warning('Empty mask at sample %d. Skipping.', m);
        continue;
    end

    %% ===================== RANDOM PHYSICAL PARAMETERS =====================
    % Intensity per valid pixel.
    indI = randi([1, NI], NL, 1);
    IA = I_max(indI);
    IA = IA(:);

    % Fiber orientation angle: 0 to 180 degrees.
    Fiber_Angle = randi([0, 360], NL, 1) ./ 2;
    Fiber_Angle = Fiber_Angle(:);

    % Generic P-SHG anisotropy parameters.
    idxP = randi([1, numel(P)], NL, 1);
    idxQ = randi([1, numel(Q)], NL, 1);

    p = P(idxP);
    q = Q(idxQ);
    p = p(:);
    q = q(:);

    % Convert p/q pair to equivalent peptide-pitch angle using lookup table.
    indPQ = sub2ind([numel(P), numel(Q)], idxP, idxQ);
    Pitch_Angle = Pitch_Angle_Table(indPQ);
    Pitch_Angle = Pitch_Angle(:);

    %% ===================== GENERIC P-SHG MODEL =====================
    I_GS = zeros(NL, N_theta);
    for mm = 1:N_theta
        I_GS(:, mm) = Generic_SHG_Model(theta(mm), Fiber_Angle, p, q);
    end

    Gain = IA ./ max(I_GS, [], 2);

    I_SHG = zeros(NX, NY, N_theta);
    for mm = 1:N_theta
        SP = I_GS(:, mm) .* Gain;
        I_SP = zeros(NX, NY);
        I_SP(Mask > 0) = SP;
        I_SHG(:, :, mm) = I_SP;
    end

    % Add Poisson noise and dark noise.
    I_SHG = fix(1e12 .* imnoise(I_SHG ./ 1e12, 'poisson')) + ...
        randi([0, Dark_Counts], NX, NY, N_theta);

    SHG_IM = single(I_SHG);

    Angle_IM = zeros(NX, NY, 2);
    FA = zeros(NX, NY);
    PA = zeros(NX, NY);

    FA(Mask > 0) = Fiber_Angle;
    PA(Mask > 0) = Pitch_Angle;

    Angle_IM(:, :, 1) = FA;
    Angle_IM(:, :, 2) = PA;
    Angle_IM = single(Angle_IM);

    save(fullfile(Save_fold, polarizationFolder, ['SHG_IM_' num2str(m) '.mat']), ...
        'SHG_IM', '-v7.3', '-nocompression');

    save(fullfile(Save_fold, angleFolder, ['Angle_IM_' num2str(m) '.mat']), ...
        'Angle_IM', '-v7.3', '-nocompression');

    clear p q Fiber_Angle Pitch_Angle I_GS I_SHG Gain SHG_IM Angle_IM FA PA;
end

fprintf('\nDone. Generic synthetic P-SHG data saved to:\n%s\n', Save_fold);
fprintf('Expected output for N = %d:\n', N);
fprintf('  Polarization_Map: %d SHG_IM_*.mat files\n', N);
fprintf('  Angle_Map       : %d Angle_IM_*.mat files\n', N);

%% ===================== LOCAL HELPER =====================
function makeFolderIfMissing(folderPath)
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end
end
