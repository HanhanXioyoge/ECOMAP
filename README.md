# ECOMAP
**Enzyme-Constrained Metabolic Modeling & Multi-Omics Analysis Platform**

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

> A comprehensive MATLAB-based toolbox for building, calibrating, and analyzing enzyme-constrained metabolic models (ecModels)

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
├── LICENSE.md                   # MIT License
│
├── images/                      # Framework diagrams and figures
│
├── scripts/                     # Core functionality modules
│   ├── Reconstruction/          # Model building and loading
│   │   ├── loadModel.m          # Unified model loader (XML/JSON/YAML/GECKO/ECMpy)
│   │   ├── convertecModel.m     # GEM to ecModel conversion
│   │   ├── fillEnzymeInformation.m
│   │   ├── getkcatfromDatabase.m
│   │   └── buildEnzConstrRxnSet.m
│   │
│   ├── Calibration/             # Model calibration
│   │   ├── GAUKS.m              # Genetic Algorithm calibration
│   │   ├── Bayesian/            # Bayesian calibration methods
│   │   │   ├── abc_max.m        # Approximate Bayesian Computation
│   │   │   ├── bayesianTuning.m
│   │   │   ├── load13CData.m
│   │   │   └── evaluateKcatRMSE.m
│   │   └── PRESTO/              # PRESTO calibration
│   │       ├── PRESTO.m
│   │       └── BatchModelgeneration.m
│   │
│   ├── Analysis/                # Model analysis and visualization
│   │   ├── ecFVA.m              # Enzyme-constrained FVA
│   │   ├── plotMultiECDF_Log10.m
│   │   └── knockout_heatmap.m
│   │
│   ├── GUI/                     # Graphical user interface
│   │   └── Core/ReconstructionApp.m
│   │
│   ├── utilities/               # Helper utilities
│   │   ├── KcatRepo.m           # kcat repository management
│   │   └── buildStoichiometricMatrix.m
│   │
│   └── ParameterManagement/     # Parameter management
│
├── tutorial/                    # MATLAB Live Script tutorials
│   ├── ecYeast_calibration_tutorial.mlx
│   ├── ecYeast_analysis_tutorial.mlx
│   ├── eciML1515_calibration_tutorial.mlx
│   ├── eciML1515_analysis_tutorial.mlx
│   ├── eciCW773_calibration_tutorial.mlx
│   └── eciCW773_analysis_tutorial.mlx
│
├── ecYeast/                    # S. cerevisiae (yeast) project
│   ├── models/                  # ecYeast9_*.mat, ecYeastGEM.yml
│   ├── data/                    # Proteomics, growth rates, flux data
│   ├── analysis/                # Bayesian states, benchmarks
│   └── ecYeast9PRESTOConfiguration.m
│
├── eciML1515/                  # E. coli iML1515 project
│   ├── models/                  # eciML1515_*.mat, iML1515.xml
│   ├── data/                    # 13C flux, proteomics, growth rates
│   ├── analysis/                # Bayesian states, benchmarks
│   └── eciML1515PRESTOConfiguration.m
│
├── eciCW773/                   # C. glutamicum project
│   ├── models/
│   ├── data/
│   └── analysis/
│
├── ecHuman/                    # Human cell model project
│
└── ecModelGEM/                 # Generic ecModel GEM framework
```

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

## Tutorial Files

The `tutorial/` directory contains MATLAB Live Scripts (.mlx) for step-by-step guidance:

| Tutorial | Description |
|----------|-------------|
| **ecYeast_calibration_tutorial.mlx** | Calibrate ecYeast model using proteomics and growth data |
| **ecYeast_analysis_tutorial.mlx** | Analyze ecYeast model predictions and benchmarks |
| **eciML1515_calibration_tutorial.mlx** | Calibrate E. coli iML1515 with multi-omics data |
| **eciML1515_analysis_tutorial.mlx** | Analyze E. coli model performance |
| **eciCW773_calibration_tutorial.mlx** | Calibrate C. glutamicum model |
| **eciCW773_analysis_tutorial.mlx** | Analyze C. glutamicum predictions |

### Running Tutorials

1. Open MATLAB
2. Navigate to the `tutorial/` folder
3. Double-click the desired `.mlx` file
4. Follow the inline instructions and run each cell sequentially

---

## Data Preparation

### Required Data Files

Each organism project (e.g., `ecYeast/`, `eciML1515/`) requires data in the `data/` subdirectory:

#### 1. Growth Rates Data

**File**: `growth_rates.tsv` or `BayesianGrowthRates.tsv`

| Column | Description | Example |
|--------|-------------|---------|
| condition | Experimental condition | Glucose, Acetate |
| growth_rate | Specific growth rate (h⁻¹) | 0.41 |
| std | Standard deviation | 0.02 |

**Example**:
```
condition	growth_rate	std
Glucose	0.41	0.02
Acetate	0.21	0.01
```

#### 2. Proteomics Data

**File**: `abs_proteomics.tsv`

| Column | Description | Example |
|--------|-------------|---------|
| gene_id | Gene identifier | b0001 |
| protein_id | UniProt ID | P0A8V2 |
| abundance | Protein abundance (fps) | 1234.5 |

#### 3. 13C Flux Data

**File**: `13CFluxdata.tsv` or `13CFluxdata.xlsx`

| Column | Description | Example |
|--------|-------------|---------|
| reaction_id | Reaction identifier | EX_glc_e |
| flux | Flux value | 10.5 |
| flux_std | Standard deviation | 0.8 |

#### 4. kcat Database (Optional)

**Directory**: `kcatData/`

**File**: `kcat_values.tsv`

| Column | Description | Example |
|--------|-------------|---------|
| ec_number | EC number | 2.7.1.1 |
| kcat | turnover number (s⁻¹) | 120 |
| source | Database source | BRENDA |

#### 5. UniProt Mapping

**File**: `uniprot.tsv`

| Column | Description | Example |
|--------|-------------|---------|
| gene_id | Gene identifier | b0001 |
| uniprot_id | UniProt ID | P0A8V2 |
| sequence | Protein sequence | MSDK... |

### Model File Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| SBML | `.xml` | Standard exchange format |
| JSON | `.json` | GECKO, ECMpy style |
| YAML | `.yml`, `.yaml` | Yeast9, ecYeastGEM |
| MATLAB | `.mat` | Various ecModel types |

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

## Project-Specific Usage

### E. coli (eciML1515)

```matlab
% Initialize project
InitializeECOMAPproject('eciML1515');

% Run full pipeline
runReconstruction;
runCalibration;
runAnalysis;
```

### Yeast (ecYeast)

```matlab
% Initialize project
InitializeECOMAPproject('ecYeast');

% Calibrate with proteomics
bayesianState = bayesianTuning('ecYeast');
```

### C. glutamicum (eciCW773)

```matlab
% Initialize project
InitializeECOMAPproject('eciCW773');
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
| `GAUKS.m` | Genetic Algorithm calibration |
| `abc_max.m` | Approximate Bayesian Computation |
| `bayesianTuning.m` | MCMC-based Bayesian tuning |
| `PRESTO.m` | Protein Reinforcement calibration |
| `load13CData.m` | Load 13C flux experimental data |

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
