# Computational Study on Baseline Arrhythmia Risk in Human Stem Cell-Derived Cardiomyocyte Phenotypes

This repository contains the anonymized simulation code used in the computational study titled "Assessing Baseline Arrhythmia Risk in Human Stem Cell-Derived Cardiomyocyte Phenotypes for Therapeutic Use in Myocardial Infarction."

## Overview

The provided code implements the validated Paci2020 electrophysiological model to simulate human induced pluripotent stem cell-derived cardiomyocytes (hPSC-CMs). The code specifically supports population-of-models (PoM) analyses for ventricular, atrial, and nodal cell phenotypes, and calculates key biomarkers relevant to arrhythmia risk assessment.

## Repository Structure

- **cellular_models/**: Contains implementations and related scripts for cellular-level simulations.
  - **Paci2020 model**: **IMPORTANT** – The original implementation of the Paci2020 model for human induced pluripotent stem cell-derived cardiomyocytes (hPSC-CMs) **is not included**. Please download it separately from the official source:
    - [Paci2020 Model Download Page](https://www.mcbeng.it/en/downloads/software/paci2020.html)

Ensure you place the downloaded Paci2020 model files within the `cellular_models/` directory before running simulations.

## Usage Instructions

1. **Download the Paci2020 model** from the provided link above.
2. **Place the model files** in the `cellular_model/` folder.
3. Execute simulations using provided scripts according to your analysis needs.

## Requirements

- MATLAB (tested on MATLAB 2020b or newer)
- No additional external toolboxes are required beyond the standard MATLAB installation.

## Disclaimer

This repository is anonymized and contains only the simulation code developed specifically for the computational study mentioned above. No experimental data or personal identifiers are included.
