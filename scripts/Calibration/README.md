# ECOMAP Calibration Module

This module calibrates enzyme-constrained models against growth, uptake and secretion fluxes, absolute proteomics measurements, and optional 13C flux data. The current implementation includes local and multi-condition sensitivity tuning, Bayesian calibration, GAUKS, PRESTO, and the sluice-parameter representation.

Return to the [project overview](../../README.md).

## Methods

| Method | Entry point | Objective | Principal requirement |
| --- | --- | --- | --- |
| Single-condition sensitivity | `sensitivityTuning` | Iteratively increase growth-limiting kcat values until a target growth rate is approached | `solveLP` |
| Multi-condition sensitivity | `MulticonditionsensitivityTuning` | Sequentially calibrate multiple carbon-source and oxygen conditions and record each kcat change | `UnconstrainedMaxGrowth.tsv` |
| Bayesian calibration | `Bayesian/bayesianTuning` | Evaluate sampled kcat sets against constrained growth, unconstrained growth, and optional 13C data | Parallel evaluation and state files |
| GAUKS | `GAUKS` | Calibrate exchange-reaction `Umin` values and protein costs on a preconstructed sluice representation | GEM, ecModel, and `UnconstrainedMaxGrowth.tsv` |
| PRESTO | `PRESTO/PRESTO` and project launch functions | Correct kcat values using absolute proteomics and growth measurements | COBRA, Gurobi, and proteomics data |
| Evaluation | `calculateCalibrationRMSE` | Report constrained-growth, unconstrained-growth, 13C, and aggregate RMSE values | The Bayesian data directory |

`applySluiceStructure`, `extractSluiceParams`, `setSluiceParams`, and `setUmin` construct and manage the sluice representation. `GAUKS` assumes that this structure has already been applied.

## Environment and project parameters

```matlab
cd('C:\path\to\ECOMAP')
setup
```

RAVEN Toolbox, COBRA Toolbox, and a functional LP/QP solver are required. The PRESTO preparation code selects Gurobi by default. Bayesian and multi-condition procedures may use Parallel Computing Toolbox.

Load project parameters as follows:

```matlab
root = fileparts(which('setup'));
manager = fullfile(root, 'projects', 'eciML1515', ...
    'eciML1515ParameterManagement.m');
params = ParameterManager.getParams(manager);

params.projectDir     = fullfile(root, 'projects', 'eciML1515');
params.calibrationDir = fullfile(params.projectDir, 'calibration');
params.analysisDir    = fullfile(params.projectDir, 'analysis');

% These fields are required by several calibration functions but are not
% currently defined by the default parameter template.
params.dataDir        = fullfile(params.calibrationDir, 'data');
params.outputDir      = params.analysisDir;
```

`bayesianTuning`, `calculateCalibrationRMSE`, `adjustKcatForAbsorptionReactions`, and PRESTO helper functions read `params.dataDir`. `GAUKS` and `MulticonditionsensitivityTuning` currently read `params.calibrationDir/data/UnconstrainedMaxGrowth.tsv` directly. The directory convention above therefore provides a consistent arrangement for all methods.

## Data organization and formats

```text
projects/<project>/calibration/
├── data/
│   ├── BayesianGrowthRates.tsv
│   ├── UnconstrainedMaxGrowth.tsv
│   ├── 13CFluxdata.tsv
│   ├── abs_proteomics.tsv
│   ├── growth_rates.tsv
│   ├── csource.tsv
│   ├── total_protein.tsv
│   └── paxDB.tsv
├── growthdata/                 # Project-specific source measurements
└── PRESTO_Results/             # PRESTO outputs
```

No single calibration method requires every file in this list.

### `BayesianGrowthRates.tsv`

Each row describes one condition, with the first column used as a row identifier. The current project data use the following logical organization:

```text
Substrate, Uptake, biomass, measured exchange reactions, OxAvail, Media
```

`Substrate` contains an exchange-reaction identifier. `OxAvail` uses `aerobic`, `anaerobic`, or `limited`. Reaction columns must match identifiers in the model; unavailable measurements should be represented by `NaN`.

### `UnconstrainedMaxGrowth.tsv`

This table is organized similarly to `BayesianGrowthRates.tsv`. The present multi-condition and GAUKS implementations also access columns 1, 3, 9, and 11 by position. Do not reorder the columns in the supplied project template: column 1 is the substrate exchange reaction, column 3 is the target growth rate, and column 11 is the oxygen condition.

### `13CFluxdata.tsv`

The recommended format is:

```text
Type  RxnName  Carbon  Direction  Cond1  Cond2 ...
```

- Rows with `Type=constraint` define biomass objectives and exchange bounds.
- Rows with `Type=flux` define 13C flux measurements. `RxnName` may contain one reaction or a grouped expression such as `{r1;r2}`.
- `Direction` uses expressions such as `{1;-1}` for grouped reactions.
- Condition-column names must begin with `Cond`.

The loader retains a legacy-format fallback, but new projects should use the explicit `Type` format.

### PRESTO inputs

| File | Organization in the current projects |
| --- | --- |
| `abs_proteomics.tsv` | UniProt or protein identifiers by row, experimental conditions by column, absolute abundance as values |
| `growth_rates.tsv` | Two columns containing condition identifiers and experimental growth rates |
| `csource.tsv` | Exchange reactions by row, conditions by column, uptake or secretion fluxes as values |
| `total_protein.tsv` | Two columns containing condition identifiers and total protein content |
| `paxDB.tsv` | Optional reference protein-abundance data |

Condition identifiers and protein identifiers must be consistent across the PRESTO inputs.

## Recommended execution sequence

### 1. Prepare the model

The model should contain:

- `enzymeConstraints.kcat`
- `enzymeConstraints.rxns`
- consistent `S`, `lb`, `ub`, and `c` fields
- valid `bioRxn` and `c_source` parameter values

After changing kcat values, call `UpdateSmatrix` to synchronize the enzyme-constraint coefficients.

### 2. Optional sensitivity calibration

```matlab
[modelTuned, changes] = sensitivityTuning( ...
    ecModel, 0.74, 5, {}, true, params);
```

The multi-condition entry point is:

```matlab
[modelTuned, summary, changes] = ...
    MulticonditionsensitivityTuning(ecModel, 5, params);
```

The multi-condition procedure returns and, when configured, writes:

- `MultiCondition_Summary.tsv`
- `MultiCondition_KcatChanges.tsv`

### 3. Bayesian calibration

```matlab
[modelBayes, finalKcats, initialKcats, rmseHistory, kcatHistory] = ...
    bayesianTuning(ecModel, true, true, true, ...
        true, 150, params.PRESTO.ncpu, 160, 0.2, [], params);
```

The three logical flags select constrained-growth data, unconstrained-growth data, and 13C data, respectively. State files are written to `params.outputDir` as `bayesianState_ABC_*.mat`; a subsequent run attempts to resume from the matching file. State files should not be reused after changing the input model, data, or sampling configuration.

### 4. GAUKS

Construct the sluice representation before calibration:

```matlab
exchangeRxns = {'EX_glc__D_e', 'EX_malt_e'};
[modelSluice, sluiceConfig] = ...
    applySluiceStructure(ecModel, exchangeRxns, 'prot_pool');

[modelAerobic, modelAnaerobic, summary] = GAUKS( ...
    modelSluice, gem, params.bioRxn, true, params);
```

The five GAUKS arguments are the ecModel, conventional GEM, GEM biomass reaction, a logical flag enabling anaerobic conditions, and the parameter structure. Earlier three-argument examples are incompatible with the current signature.

### 5. PRESTO

The low-level `PRESTO` function receives a model cell array, experimental growth rates, a protein-abundance matrix, and optional correction parameters. Project workflows should preferentially use `startPRESTOproject` or `startPRESTO` to prepare these inputs from the project configuration. Before execution, confirm that:

- `params.PRESTO.ncpu`, `nIter`, `epsilon`, `lambda`, and `theta` are defined;
- condition identifiers agree across all four TSV inputs;
- COBRA is configured with a functional Gurobi LP solver; and
- the protein prefixes in the ecModel agree with the PRESTO configuration.

## Web bridges

| Operation | Bridge | Optimization class |
| --- | --- | --- |
| Apply sluice structure | `mdpApplySluice` | None |
| Initialize kcat repository | `mdpKcatRepoInit` | None |
| Sensitivity tuning | `mdpSensitivityTuning` | Underlying LP operations |
| GAUKS | `mdpGauks` | LP |
| Bayesian calibration | `mdpBayesian` | LP and parallel sample evaluation |
| PRESTO | `mdpPresto` | QP/LP |

Web jobs map MATLAB failures to the common bridge catalogue, including `err_no_proteomics`, `err_sluice_data`, `err_presto_data`, and `err_gurobi_license`.

## Troubleshooting

### `Reference to non-existent field 'dataDir'` or `outputDir`

Define the aliases shown under “Environment and project parameters.” The current `Template.m` defines stage directories but does not define these two fields.

### Target, biomass, or carbon-source reactions cannot be found

`params.bioRxn`, `params.c_source`, the `Substrate` entries, and all flux-column names must use reaction identifiers present in the model.

### Multi-condition results appear misaligned

Retain the column ordering of the supplied `UnconstrainedMaxGrowth.tsv` template. The current multi-condition and GAUKS implementations use both variable names and fixed column indices.

### Bayesian calibration resumes an obsolete run

Inspect and relocate the relevant `bayesianState_ABC_*.mat` file. These files contain model-specific sampling state and should not be shared across models or datasets.

### PRESTO or Gurobi fails

Verify the Gurobi MATLAB interface, license, and COBRA solver configuration. Alternative LP solvers may support other ECOMAP procedures, but the current PRESTO preparation explicitly selects `gurobi`.
