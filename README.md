# ECOMAP
**Enzyme-Constrained Optimization of Metabolic Models and Analysis Pipeline**

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

> A comprehensive MATLAB-based toolbox for reconstructing, calibrating, and analyzing enzyme-constrained metabolic models (ecModels)

## Framework Overview

![ECOMAP Framework](images/ECOMAP_workflow.png)

ECOMAP provides a unified framework for:

1. **Model Reconstruction** - Convert traditional GEMs to enzyme-constrained models by integrating kcat values and protein constraints
2. **Model Calibration** - Refine model parameters using proteomics, growth rates, and 13C flux data
3. **Model Analysis** 

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
└── ecHuman/                    # Human cell model project
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
| **growth_rates.tsv** | Maximum growth rate constraints | `PRESTO` | Tab-separated: condition name, max growth rate |
| **BayesianGrowthRates.tsv** | Growth rates with detailed flux data | `bayesianTuning` (useConstraint=true) | Tab-separated with substrate, uptake, and reaction fluxes |
| **UnconstrainedMaxGrowth.tsv** | Growth rates without substrate constraints | `bayesianTuning` (useUnconstrained=true), `GAUKS` | Tab-separated with substrate and reaction fluxes |
| **13CFluxdata.tsv** | 13C metabolic flux data | `bayesianTuning` (use13Cflux=true) | Tab-separated: reaction ID, flux value, std |
| **csource.tsv** | Carbon source exchange reaction fluxes | `PRESTO` | Matrix format: reactions × conditions |

### Multi-Omics Data

| File | Purpose | Calibration Method | Format |
|------|---------|-------------------|--------|
| **abs_proteomics.tsv** | Absolute protein abundances (fps) for PRESTO | `PRESTO` | Gene ID, UniProt ID, abundance |
| **paxDB.tsv** | Protein abundance database for cross-referencing | Model annotation | Gene ID, abundance score |
| **total_protein.tsv** | Total protein content for protein pool calculation | `updateProtPool` | Total protein (g/gDCW) |

### Model Annotation Data

| File | Purpose | Used By | Format |
|------|---------|---------|--------|
| **uniprot.tsv** | UniProt ID mapping and protein sequences | `fillEnzymeInformation`, `writeInputFile` | Gene ID, UniProt ID, sequence |
| **metInfo.tsv** | Metabolite information (SMILES, InChI, etc.) | `getMetinfo` | Metabolite ID, identifiers |
| **ComplexPortal.json** | Protein complex annotations | `applyComplexdata` | JSON format |
| **kcatData/** (folder) | kcat values from databases (BRENDA, SABIO) | `getkcatfromDatabase`, `dbKcatMatch` | TSV with EC number, kcat, source |

### Data File Formats

#### growth_rates.tsv (for PRESTO)
Simple two-column format: condition name and maximum growth rate.
```
ConditionName    MaxGrowthRate
DiBartolomeo2020_Gluc    0.398366667
DiBartolomeo2020_Etoh    0.1225
Lahtvee2017_REF    0.1
Yu2021_N30_005    0.05
```

#### BayesianGrowthRates.tsv (for bayesianTuning with useConstraint=true)
Multi-column format with substrate uptake and reaction fluxes.
```
    Substrate    Uptake    r_2111    r_1634    r_1761    r_1808    r_2033    r_1672    r_1992    r_1654    OxAvail    Media
1    r_1714    -13.33333333    0.36    NaN    19.82608696    0.38    0    NaN    NaN    NaN    aerobic    MIN
2    r_1714    -15.44444444    0.39    NaN    22.4068323    0.37    0    NaN    NaN    NaN    aerobic    MIN
...
```
- **Substrate**: Reaction ID for carbon source
- **Uptake**: Substrate uptake rate (negative = consumption)
- **r_XXXX**: Flux through various reactions (r_2111 = growth rate)
- **OxAvail**: aerobic/anaerobic/limited
- **Media**: MIN or other media types

#### UnconstrainedMaxGrowth.tsv (for bayesianTuning with useUnconstrained=true)
Similar to BayesianGrowthRates but only contains growth rates (no uptake data).
```
    Substrate    Uptake    r_2111    r_1634    r_1761    r_1808    r_2033    r_1672    r_1992    r_1654    OxAvail    Media
1    r_1709    NaN    0.338    NaN    NaN    NaN    NaN    NaN    NaN    NaN    aerobic    MIN
2    r_1710    NaN    0.28    NaN    NaN    NaN    NaN    NaN    NaN    NaN    aerobic    MIN
3    r_1714    NaN    0.41    NaN    NaN    NaN    NaN    NaN    NaN    NaN    aerobic    MIN
...
```

#### csource.tsv (for PRESTO getconditions)
Matrix format with exchange reactions as rows and experimental conditions as columns.
```
exchangeRxn    DiBartolomeo2020_Gluc    DiBartolomeo2020_Etoh    Lahtvee2017_REF    ...
r_1714    -14.38395556    0    -1.070928479    ...
r_1808    0.666410217    -0.202562355    0.000857154    ...
r_1634    0.444200167    0.0515425    0    ...
r_1992    -1000    -1000    -2.446472988    ...
r_1672    0    0    2.62755453    ...
r_1654    -1000    -1000    -1000    ...
...
```

#### abs_proteomics.tsv
```
gene_id    protein_id    abundance
YAL038W    P00560    1234.5
YBR019C    P32167    567.8
```

#### 13CFluxdata.tsv
13C metabolic flux data with constraint and flux types.
```
Type       RxnName                  Carbon    Direction    Cond1     Cond2     Cond3     Cond4
constraint r_2111                                0.405      0.150     0.300     0.400
constraint r_1714                               -16.731     -1.560     -4.900     -8.230
constraint r_1672                                28.759      3.610     10.320    13.740
constraint r_1992                                -3.269      Nan       Nan       Nan
flux      r_0534                    6          {1}        16.731     1.560     4.900     8.230
flux      r_0466                    6          {1}         1.781     0.853     1.450     1.243
flux      r_0467                    6          {1}        14.168     0.393     2.940     6.411
...
```
- **Type**: `constraint` (growth/substrate constraints) or `flux` (reaction flux values)
- **RxnName**: Reaction identifier
- **Carbon**: Number of carbon atoms (only for flux type)
- **Direction**: Reaction directionality notation (e.g., {1}, {-1;1})
- **Cond1-Cond4**: Flux values under different experimental conditions

#### kcatData/kcat_values.tsv
```
ec_number    kcat    source
2.7.1.1    120    BRENDA
1.1.1.1    85    SABIO
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