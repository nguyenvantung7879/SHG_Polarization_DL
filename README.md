# Physics-pretrained transfer learning for P-SHG collagen mapping

MATLAB code, trained networks, experimental examples, and saved results for mapping **collagen fiber orientation** and **peptide-pitch angle** from polarization-resolved second-harmonic generation (P-SHG) microscopy.

The workflow combines physics-based synthetic pretraining with fine-tuning on 64 experimental patches. An 18-state P-SHG stack is the input; two spatial angle maps, in degrees, are the outputs.

## System overview

![Figure 1: synthetic pretraining and experimental transfer-learning workflow](figure/Figure%201.png)

*Figure 1. Synthetic P-SHG generation from digit masks and forward SHG models, followed by transfer learning using experimental polarization stacks and pixel-wise fitting labels. The diagram also shows the dual-liquid-crystal P-SHG acquisition system.*

```text
Synthetic P-SHG stacks + angle labels -> synthetic-pretrained network
64 experimental input/label patches + pretrained network -> fine-tuned network
18-state P-SHG stack -> network -> fiber orientation + peptide pitch
```

## Representative tissue results

![Figure 6: intensity, fiber-orientation and peptide-pitch maps, and angle distributions](figure/Figure%206.png)

*Figure 6. Representative subchondral bone (top), native cartilage (middle), and regenerative tissue (bottom). Each row shows P-SHG intensity, fiber orientation and its distribution, and peptide-pitch angle and its distribution. Click the image to inspect the full-resolution figure.*

This is the archived manuscript figure. Its regenerative inputs were clipped per polarization frame at the 90th percentile and rescaled by 13; the bone and cartilage inputs were not processed this way. Interpret the comparison with this preprocessing difference in mind. Histogram peaks are descriptive summaries, not pixel-wise accuracy measurements.

## Contents

```text
figure/                        Figure 1 and Figure 6 manuscript images
model/                         Three trained MATLAB network checkpoints
function/                      Inference, forward-model, and plotting helpers
data_test/
  Polarizer/                   Multipage experimental P-SHG TIFF inputs
  Angle/                       Pixel-wise fitting reference TIFFs
  Tendon_Model_Evaluation_Results/
  Tissue_Model_Evaluation_Results/
data_transfer_learning/
  Selected_64_data/
    Polarization_32x32_Mat/     64 input patches (32 x 32 x 18)
    Angle_32x32_Mat/            64 paired labels (32 x 32 x 2)
    selected_patch_list.csv    Patch-selection/source metadata
data_simulation/               Synthetic generator, forward models, lookup files
run_tendon_eval.m              Compare all three networks on tendon data
run_tissue_eval.m              Proposed-model tissue evaluation and Figure 6-style plots
```

The test inputs contain 20 tendon stacks and 10 stacks each for `Col I`, `Col II`, and `Regerative tissue`. Keep the existing spelling `Regerative tissue`, which is used by the scripts.

## Requirements and installation

The original workstation used MATLAB R2025b on Windows, an Intel Core i5-12400, 32 GB RAM, and an NVIDIA GeForce RTX 3060 with 12 GB memory.

| Requirement | Use |
|---|---|
| MATLAB R2025b | Development environment and TIFF/MAT-file processing |
| Deep Learning Toolbox | Load networks and compute activations |
| Image Processing Toolbox | Image processing, overlays, synthetic generation |
| Curve Fitting Toolbox | Histogram-peak fitting and lookup generation |
| Parallel Computing Toolbox and a supported NVIDIA GPU | GPU inference in the supplied helper |

```powershell
git clone https://github.com/nguyenvantung7879/SHG_Polarization_DL.git
cd SHG_Polarization_DL
```

Open MATLAB and set its **Current Folder to this repository root**, containing `run_tendon_eval.m` and `run_tissue_eval.m`. The evaluation scripts use paths relative to this folder and Windows path separators.

```matlab
addpath(fullfile(pwd, 'function'));
assert(canUseGPU, 'The supplied inference helper requires a supported GPU.');
gpuDevice
```

No training or synthetic generation is required to use the bundled checkpoints.

## Trained models

Each MAT checkpoint contains a network variable named `net`.

| Model | Checkpoint | Role |
|---|---|---|
| Synthetic-pretrained only | [SHG_Virtual_MNIST_Digit-IMAGE_Generic_Angle-ResNet-Filter-128-Epoch-200-Time-284624.mat](model/SHG_Virtual_MNIST_Digit-IMAGE_Generic_Angle-ResNet-Filter-128-Epoch-200-Time-284624.mat) | Synthetic pretraining without experimental fine-tuning |
| Without pretraining | [TrainWithoutPretrain_64RealPatch-Epoch-200-Selected32Patch-64-TrainPatchSize-16-Balance-ColI50-ColII50-Seed-7-Time-260.mat](model/TrainWithoutPretrain_64RealPatch-Epoch-200-Selected32Patch-64-TrainPatchSize-16-Balance-ColI50-ColII50-Seed-7-Time-260.mat) | Experimental training on 64 patches without synthetic pretraining |
| Proposed: pretrained + fine-tuned | [Finetune_Proposed_64RealPatch-Epoch-200-Selected32Patch-64-TrainPatchSize-16-Balance-ColI50-ColII50-Seed-7-Time-260.mat](model/Finetune_Proposed_64RealPatch-Epoch-200-Selected32Patch-64-TrainPatchSize-16-Balance-ColI50-ColII50-Seed-7-Time-260.mat) | Synthetic pretraining followed by experimental fine-tuning; used for tissue evaluation |

The helper reads the `Angle_Predict` layer: channel 1 is fiber orientation and channel 2 is peptide pitch. Outputs retain their checkpoint's label convention. The original synthetic pitch convention differs from the experimental fitting convention; account for that distinction when comparing models.

## Predict one stack

Run this example in MATLAB from the repository root. It uses an included tendon input and does not overwrite saved evaluation results.

```matlab
addpath(fullfile(pwd, 'function'));
modelFiles = dir(fullfile('model', 'Finetune_Proposed_64RealPatch-Epoch-200-*.mat'));
assert(numel(modelFiles) == 1, 'Expected exactly one proposed checkpoint.');
model = load(fullfile(modelFiles(1).folder, modelFiles(1).name), 'net');

inputPath = fullfile('data_test', 'Polarizer', 'Tendon', 'Tendon_1_ROI01_TL.tif');
stack = single(tiffreadVolume(inputPath));
assert(size(stack, 3) >= 18, 'Expected at least 18 polarization states.');
stack = stack(:, :, 1:18);

threshold = 150;
mask = sum(stack, 3) <= threshold;
[fiber, pitch, inferenceSeconds] = predict_image(model.net, stack, threshold);
fiber = double(gather(fiber));
pitch = double(gather(pitch));
fiber(mask) = NaN;
pitch(mask) = NaN;

figure;
tiledlayout(1, 2);
nexttile; imagesc(fiber, [0 180]); axis image off; colorbar;
title('Fiber orientation (degrees)');
nexttile; imagesc(pitch, [40 60]); axis image off; colorbar;
title('Peptide pitch (degrees)');
fprintf('Inference: %.3f seconds\n', inferenceSeconds);
```

`predict_image` always requests GPU execution. Its `Threshold` argument is unused internally, so the example explicitly masks the returned 2-D maps. `[40 60]` is a display window, not a constraint on raw pitch predictions.

## Run the evaluation workflows

Both scripts write into `data_test/*_Model_Evaluation_Results/` and may overwrite archived outputs. Copy those folders first if you need to preserve the existing results.

### Rat tendon: compare three models

```matlab
run_tendon_eval
```

This pairs the TIFFs in `data_test/Polarizer/Tendon/` with same-named references in `data_test/Angle/Tendon/`, loads the three checkpoints, and exports maps, histogram-peak/total-variation summaries, and representative plots.

Outputs in [Tendon_Model_Evaluation_Results](data_test/Tendon_Model_Evaluation_Results/):

- `tendon_eval_paper_table.csv`, `tendon_eval_numeric_summary.csv`, `tendon_eval_per_sample.csv`.
- `Prediction_TIFs/`: two-page prediction TIFFs for each model.
- `Representative_Plots/`: PNG and MATLAB FIG plots.

Settings at the top control thresholds, sample exclusions, the representative ROI, and output switches. Set `excludeSampleKeywords = {};` to explicitly include all samples.

### Porcine tissue: proposed model

```matlab
run_tissue_eval
```

This processes `Col I`, `Col II`, and `Regerative tissue`. `fileIndexBySample = [10, 10, 1]` selects the displayed stacks after alphabetical sorting; `fileKeywordBySample` can select them by name. The script also scans all available tissue TIFFs for descriptive tables.

Outputs in [Tissue_Model_Evaluation_Results](data_test/Tissue_Model_Evaluation_Results/):

- `Fig6_ColI_ColII_Reg_per_file.csv` and `Fig6_ColI_ColII_Reg_summary_by_sample.csv`.
- `Fig6_selected_figure_peaks.csv`: displayed files and histogram peaks.
- `Prediction_TIFs/` and `Representative_Plots/`.

The script preserves the archived processing: its 2-D logical indexing into the 3-D input stack affects only the first polarization plane, and descriptive-table masks differ from selected-figure masks. `ValidPixelFraction` does not measure coverage above the configured intensity threshold. Re-running generates the current script's plots; the manuscript layout is provided separately as `figure/Figure 6.png`.

## Data formats

| Data | Format |
|---|---|
| Experimental input | Multipage TIFF, `H x W x 18`, ordered `0:10:170` degrees; evaluation scripts retain the first 18 pages |
| Fitting reference | Same-named two-page TIFF: fiber orientation on page 1, peptide pitch on page 2, in degrees |
| Saved prediction | Two-page angle TIFF in the same channel order |
| Fine-tuning input / target | MAT arrays `32 x 32 x 18` / `32 x 32 x 2` |

Use the intensity/count scale appropriate to the checkpoint. Arbitrary normalization or per-frame clipping changes its input. Fitting references are model-derived estimates, not independent ground truth.

## Synthetic generation and training data

[data_simulation](data_simulation/) contains `Generate_Virtual_Polarization_Spectrum_Digit_Image.m`, `generate_pitch_lookup_pq.m`, both forward models, and the precomputed `Pitch_Angle_Lookup.mat` / `P_Q_Value.mat`. The generator creates digit-mask stacks and paired labels with photon/dark noise. Its final `N = 1` assignment selects a one-image demonstration.

These scripts retain paths from an earlier layout. Before running the generator in this checkout, update these assignments:

```matlab
repoRoot = fileparts(scriptFolder);
Save_fold = fullfile(repoRoot, 'data_simulation', 'Synthetic_PSHG_MNIST');
lookupPath = fullfile(repoRoot, 'data_simulation', 'Pitch_Angle_Lookup.mat');
pqPath = fullfile(repoRoot, 'data_simulation', 'P_Q_Value.mat');
```

Set the final `N` assignment to the intended demonstration size. The lookup files are bundled, so regeneration is optional. To run `generate_pitch_lookup_pq.m`, set `repoRoot = fileparts(scriptFolder)` and `modelFolder = fullfile(repoRoot, 'data_simulation')`. The generator also requires the MATLAB `DigitDataset` example folder checked by its code.

The 64 experimental input/target pairs and source metadata are included in [Selected_64_data](data_transfer_learning/Selected_64_data/). Complete synthetic-pretraining and transfer-training entry points are not bundled in this checkout; the workflows above use the saved checkpoints.

## Data availability and license

This repository supplies the test stacks, selected fine-tuning patches, trained checkpoints, fitting references, and saved results. It is not a complete archive of all raw imaging sessions. Additional raw data may be requested from the corresponding authors.

The source code is released under the [MIT License](LICENSE). Please cite the associated manuscript when using this work.
