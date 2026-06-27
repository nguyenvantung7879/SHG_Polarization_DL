# P-SHG Collagen Physics-Pretrained Transfer Learning

This repository contains MATLAB code and representative data for rapid spatial quantification of collagen fiber orientation and peptide-pitch angle from polarization-resolved second-harmonic generation (P-SHG) microscopy images using physics-pretrained transfer learning.

The workflow includes physics-based synthetic P-SHG data generation, transfer learning with limited experimental data, and evaluation on rat tendon and porcine knee tissue samples.

## Project overview

Conventional pixel-wise fitting of P-SHG images can estimate collagen structural parameters, but it is computationally expensive and often produces fragmented angle maps. This project implements a deep learning-based framework that directly predicts two spatial maps from an 18-state P-SHG image stack:

* collagen fiber orientation angle
* peptide-pitch angle

The model is pretrained using synthetic P-SHG data generated from a physics-based generic SHG model and then fine-tuned using a limited number of experimentally labeled P-SHG patches.

## Repository structure

```text
.
├── data_simulation/
│   └── Generated code
│      
├── data_test/
│   ├── Polarizer/
│   └── Angle/
│
├── data_transfer_learning/
│   └── Selected_64_data/
│
├── function/
│   ├── Generic_SHG_Model.m
│   ├── predict_image.m
│   └── other helper functions
│
├── model/
│   ├── Pitch_Angle_Lookup.mat
│   ├── P_Q_Value.mat
│   └── trained model files
│
├── run_tendon_eval.m
├── run_tissue_eval.m
├── LICENSE
└── README.md
```

## Software requirements

The code was developed and tested using:

```text
MATLAB R2025b
Windows 11 64-bit
Deep Learning Toolbox
Image Processing Toolbox
Curve Fitting Toolbox
```

The workstation used for model training and inference was:

```text
CPU: Intel Core i5-12400
RAM: 32 GB
GPU: NVIDIA GeForce RTX 3060, 12 GB
```

## Notes on data availability

This repository includes analysis scripts, model evaluation code, processed result tables, and representative data for reproducibility.

Full raw P-SHG imaging datasets may be large and are not necessarily included in this repository. Additional raw imaging data are available from the corresponding authors upon reasonable request.


## License

The source code in this repository is released under the MIT License.

Representative data and processed results are provided for academic and reproducibility purposes. Please cite the associated manuscript if you use this repository.
