# ECOMAP
**Enzyme-Constrained Metabolic Modeling & Multi-Omics Analysis Platform**

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

> A comprehensive MATLAB-based toolbox for reconstructing, calibrating, and analyzing enzyme-constrained metabolic models (ecModels)

## Framework Overview

![ECOMAP Framework](images/ECOMAP_workflow.png)

ECOMAP provides a unified framework for:

1. **Model Reconstruction** - Convert traditional GEMs to enzyme-constrained models by integrating kcat values and protein constraints
2. **Model Calibration** - Refine model parameters using proteomics, growth rates, and 13C flux data
3. **Model Analysis** - Evaluate predictions through FVA, gene knockout studies, and protein efficiency assessment

---

## Directory Structure

```
ECOMAP/
├── README.md                    # This file
├── setup.m                      # MATLAB path setup script
├── images/                      # Framework diagrams
│
├── scripts/                     # Core functionality
│   ├── Reconstruction/          # Model building and loading
│   ├── Calibration/             # Bayesian, GAUKS, PRESTO methods
│   ├── Analysis/                # FVA, visualization, knockout
│   ├── GUI/                     # Graphical user interface
│   ├── utilities/               # KcatRepo, model I/O
│   └── ParameterManagement/     # Parameter management system
│
├── tutorial/                    # MATLAB Live Script tutorials (see below)
│
├── ecYeast/                    # S. cerevisiae project
├── eciML1515/                  # E. coli iML1515 project
├── eciCW773/                   # C. glutamicum project
├── ecHuman/                    # Human cell model project
└── ecModelGEM/                 # Generic ecModel GEM framework
```

Each organism project follows the structure:
```
[project]/
├── models/           # Model files (.mat, .xml, .yml)
├── data/             # Experimental data (see Data Preparation)
└── analysis/         # Calibration results and benchmarks
```

---

## Tutorial Files

The `tutorial/` directory contains MATLAB Live Scripts (.mlx) for step-by-step guidance:

| Tutorial | Description |
|----------|-------------|
| **tutorial.mlx** | General introduction to ECOMAP |
| **ecYeast_reconstruction_tutorial.mlx** | Build ecYeast model from GEM |
| **ecYeast_calibration_tutorial.mlx** | Calibrate ecYeast with proteomics and growth data |
| **ecYeast_analysis_tutorial.mlx** | Analyze ecYeast predictions and benchmarks |
| **eciML1515_reconstruction_tutorial.mlx** | Build E. coli iML1515 ecModel |
| **eciML1515_calibration_tutorial.mlx** | Calibrate E. coli with Bayesian/PRESTO methods |
| **eciML1515_analysis_tutorial.mlx** | Analyze E. coli model performance |
| **eciCW773_reconstruction_tutorial.mlx** | Build C. glutamicum ecModel |
| **eciCW773_calibration_tutorial.mlx** | Calibrate C. glutamicum model |
| **eciCW773_analysis_tutorial.mlx** | Analyze C. glutamicum predictions |
| **ecHuman_tutorial.mlx** | Human cell model reconstruction and analysis |

### Running Tutorials

1. Open MATLAB
2. Navigate to the `tutorial/` folder
3. Double-click the desired `.mlx` file
4. Follow the inline instructions and run each cell sequentially

---

## Data Preparation

Each organism project (e.g., `ecYeast/`, `eciML1515/`) requires data files in the `data/` subdirectory. Below is the complete list based on the ecYeast project structure.

### Core Calibration Data

| File | Purpose | Calibration Method | Format |
|------|---------|-------------------|--------|
| **growth_rates.tsv** | Growth rate measurements for sensitivity analysis | `sensitivityTuning`, `MulticonditionsensitivityTuning` | Tab-separated with condition name and growth rate |
| **BayesianGrowthRates.tsv** | Growth rate data for Bayesian parameter estimation | `bayesianTuning` (ABC method) | Tab-separated with substrate and growth rate |
| **UnconstrainedMaxGrowth.tsv** | Maximum growth rate constraints for GAUKS calibration | `GAUKS` | Tab-separated: condition name, max growth rate |
| **csource.tsv** | Carbon source information and uptake bounds | PRESTO configuration | Tab-separated |

### Multi-Omics Data

| File | Purpose | Calibration Method | Format |
|------|---------|-------------------|--------|
| **abs_proteomics.tsv** | Absolute protein abundances (fps) for PRESTO | `PRESTO`, `getconditions` | Gene ID, UniProt ID, abundance |
| **paxDB.tsv** | Protein abundance database for cross-referencing | Model annotation | Gene ID, abundance score |
| **total_protein.tsv** | Total protein content for protein pool calculation | `updateProtPool` | Total protein (g/gDCW) |

### Flux Data

| File | Purpose | Calibration Method | Format |
|------|---------|-------------------|--------|
| **13CFluxdata.tsv** | 13C metabolic flux data for flux balance constraints | `bayesianTuning` (use13Cflux=true) | Tab-separated: reaction ID, flux value, std |
| **growthdata/** (folder) | Substrate uptake and growth data for verification | Model validation | Contains Aerobic.tsv, Anaerobic.tsv, Mul_csources.tsv |

### Model Annotation Data

| File | Purpose | Used By | Format |
|------|---------|---------|--------|
| **uniprot.tsv** | UniProt ID mapping and protein sequences | `fillEnzymeInformation`, `writeInputFile` | Gene ID, UniProt ID, sequence |
| **metInfo.tsv** | Metabolite information (SMILES, InChI, etc.) | `getMetinfo` | Metabolite ID, identifiers |
| **ComplexPortal.json** | Protein complex annotations | `applyComplexdata` | JSON format |
| **kcatData/** (folder) | kcat values from databases (BRENDA, SABIO) | `getkcatfromDatabase`, `dbKcatMatch` | TSV with EC number, kcat, source |

### Data File Formats

#### growth_rates.tsv / BayesianGrowthRates.tsv
```
Substrate    GrowthRate
Glucose      0.41
Acetate      0.21
Ethanol      0.12
```

#### UnconstrainedMaxGrowth.tsv
```
ConditionName    MaxGrowthRate
Glucose          0.41
Ethanol          0.12
```

#### abs_proteomics.tsv
```
gene_id    protein_id    abundance
YAL038W    P00560       1234.5
YBR019C    P32167       567.8
```

#### 13CFluxdata.tsv
```
reaction_id    flux    flux_std
EX_glc_e       10.5    0.8
BIOMASS        0.41    0.02
```

#### kcatData/kcat_values.tsv
```
ec_number    kcat    source
2.7.1.1      120     BRENDA
1.1.1.1      85      SABIO
```

### Model File Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| SBML | `.xml` | Standard exchange format |
| JSON | `.json` | GECKO, ECMpy style |
| YAML | `.yml`, `.yaml` | Yeast9, ecYeastGEM |
| MATLAB | `.mat` | Various ecModel types |

---

## Installation

### Requirements

- **MATLAB** 2022a or higher
- **COBRA Toolbox** 3.1+ ([Installation guide](https://opencobra.github.io/cobratoolbox/stable.html))
- **RAVEN Toolbox** 2.9.2+ ([Installation guide](https://sysbiochalmers.github.io/RAVEN/))
- **Gurobi Optimizer** (or other MATLAB-supported solvers)

### Setup

1. Clone the repository:
```bash
git clone https://github.com/HanhanXioyoge/ECOMAP.git
cd ECOMAP
```

2. Add all paths in MATLAB:
```matlab
setup
```

3. Verify installation:
```matlab
>> which loadModel.m
>> which convertecModel.m
```

---

## Quick Start

### 1. Load a Model

```matlab
% Load from SBML
model = loadModel('iML1515.xml');

% Load from YAML
model = loadModel('ecYeastGEM.yml');

% Load from MAT
model = loadModel('eciML1515_PRESTO.mat');
```

### 2. Convert GEM to ecModel

```matlab
% Convert to enzyme-constrained model
ecModel = convertecModel(model, 'kcatSource', 'BRENDA');
```

### 3. Calibrate Model

```matlab
% Bayesian calibration
[bayesianState, RMSE] = abc_max(projectFolder);

% PRESTO calibration
[calibratedModel, params] = PRESTO(projectFolder);
```

### 4. Analyze Model

```matlab
% Flux variability analysis
[fmin, fmax] = ecFVA(ecModel);

% Multi-condition analysis
plotMultiECDF_Log10(resultsFolder);
```

---

## Key Functions Reference

### Reconstruction
| Function | Description |
|----------|-------------|
| `loadModel.m` | Load models (XML, JSON, YAML, MAT) |
| `convertecModel.m` | Convert GEM to ecModel |
| `fillEnzymeInformation.m` | Add UniProt IDs and sequence info |
| `getkcatfromDatabase.m` | Retrieve kcat from BRENDA/SABIO |
| `buildEnzConstrRxnSet.m` | Build enzyme constraint reactions |

### Calibration
| Function | Description |
|----------|-------------|
| `GAUKS.m` | Growth-Anchored kcat calibration |
| `abc_max.m` | Approximate Bayesian Computation |
| `bayesianTuning.m` | MCMC-based Bayesian tuning |
| `PRESTO.m` | Protein Reinforcement calibration |
| `sensitivityTuning.m` | Single-condition sensitivity analysis |
| `MulticonditionsensitivityTuning.m` | Multi-condition sensitivity analysis |

### Analysis
| Function | Description |
|----------|-------------|
| `ecFVA.m` | Enzyme-constrained FVA |
| `plotEcFVA.m` | Visualize FVA results |
| `knockout_heatmap.m` | Gene knockout analysis |
| `plotMultiECDF_Log10.m` | Multi-condition CDF plots |
| `plotRMSEBar.m` | RMSE comparison plots |

---

## Output Files

After calibration/analysis, results are stored in each project's `analysis/` folder:

| File | Description |
|------|-------------|
| `bayesianState_*.mat` | Bayesian posterior states |
| `MultiCondition_Summary.tsv` | Multi-condition results |
| `Benchmark_*.png` | Prediction benchmarks |
| `*_kcatRepo.mat` | kcat repository |
| `MultiCondition_KcatChanges.tsv` | kcat parameter changes |

---

## License

This project is licensed under the MIT License. See [LICENSE.md](LICENSE.md) for details.

---

## Citation

If you use ECOMAP in your research, please cite:

```
ECOMAP: Enzyme-Constrained Metabolic Modeling & Multi-Omics Analysis Platform
```