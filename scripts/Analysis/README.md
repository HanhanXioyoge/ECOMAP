# ECOMAP Analysis Module

This module provides flux-variability analysis for enzyme-constrained models, statistical visualization of kcat and RMSE results, and Live Scripts for fermentation and overflow-metabolism studies. The web analysis workflow additionally exposes gene-knockout and protein-usage aggregation through MATLAB bridges.

Return to the [project overview](../../README.md).

## Current contents

| File | Type | Purpose |
| --- | --- | --- |
| `ecFVA.m` | Function | Perform FVA on a GECKO 3-style ecModel and map the result to the reaction order of its GEM counterpart |
| `plotEcFVA.m` | Function | Plot cumulative distributions of flux-variability ranges for one or more models |
| `plotMultiECDF_Log10.m` | Function | Plot log10 empirical cumulative distributions and median guides for multiple datasets |
| `plotRMSEBar.m` | Function | Produce sortable, publication-oriented RMSE bar plots |
| `plotSciPie.m` | Function | Summarize categorical data and generate a scientific-color pie chart |
| `chemostat_batch_simulation.mlx` | Live Script | Chemostat and batch-culture simulation |
| `overflow_metabolism.mlx` | Live Script | Overflow-metabolism analysis |
| `kcat_up.mlx` | Live Script | kcat adjustment and result inspection |

Files referenced by earlier documentation—`assessProteinEfficiency.mlx`, `knockout_heatmap.m`, and `predict_plot.mlx`—are not present in this directory. Some corresponding capabilities now reside in web bridges or other scripts and cannot be invoked under those former filenames.

## Dependencies

- Run `setup` from the repository root.
- RAVEN Toolbox and a model-toolbox environment providing `solveLP`, `setParam`, and `mapRxnsToConv`.
- A functional LP solver.
- The direct `ecFVA` implementation creates a parallel pool and calls `progressbar`. It therefore requires Parallel Computing Toolbox and a resolvable `progressbar` implementation. If this parallel path fails, the web bridge attempts a serial FVA implementation.

## ecFVA

### Inputs and outputs

```matlab
[minFlux, maxFlux] = ecFVA(ecModel, gemModel);
```

- `ecModel`: a GECKO 3-style enzyme-constrained model containing reactions, stoichiometry, bounds, an objective, and the relevant `ec` or enzyme-constraint structures.
- `gemModel`: the corresponding unconstrained GEM. The final outputs are mapped to the order of `gemModel.rxns`.

The function consolidates reaction variants with `_REV` and `_EXP_n` suffixes and prevents enzyme-usage reactions from contaminating the final GEM-level mapping. The returned `minFlux` and `maxFlux` vectors are aligned with `gemModel.rxns`.

### Minimal example

```matlab
setup

ecData = load('projects/eciML1515/models/eciML1515-integrated.mat');

% Variable names are not uniform across all MAT files.
fieldnames(ecData)

root = fileparts(which('setup'));
manager = fullfile(root, 'projects', 'eciML1515', ...
    'eciML1515ParameterManagement.m');
params = ParameterManager.getParams(manager);
params.modelsDir = fullfile(root, 'projects', 'eciML1515', 'models');

ecModel = ecData.ecModel;  % Replace with the actual variable name.
gemModel = loadModel('iML1515.xml', 'TRADITION', params.modelsDir, params);
[minFlux, maxFlux] = ecFVA(ecModel, gemModel);
plotEcFVA(minFlux, maxFlux);
```

The explicit `fieldnames` inspection is necessary because the repository does not enforce a single variable name across all `.mat` artifacts.

### Multi-model comparison

Provide one result per column:

```matlab
allMin = [minFluxModelA, minFluxModelB];
allMax = [maxFluxModelA, maxFluxModelB];
plotEcFVA(allMin, allMax);
```

Such a comparison is meaningful only if all columns have already been mapped to the same reaction ordering.

## Visualization utilities

### Multiple empirical cumulative distributions

```matlab
out = plotMultiECDF_Log10( ...
    {kcatBefore, kcatAfter}, ...
    {'Before', 'After'}, ...
    'Title', 'kcat distribution', ...
    'XLabel', 'kcat (s^{-1})');
```

Input may be a cell array or a numerical matrix with one dataset per column. The result structure contains the figure, axes, line objects, filtered data, and the 50th-percentile intersection for each group.

### RMSE bar plot

```matlab
out = plotRMSEBar( ...
    {'Initial', 'Bayesian', 'PRESTO'}, ...
    [0.31, 0.18, 0.12], ...
    'Sort', 'descend', ...
    'YLabel', 'RMSE');
```

### Categorical pie chart

```matlab
[fig, ax, stats] = plotSciPie(kcatSources, ...
    'Title', 'kcat source', ...
    'Order', 'descend', ...
    'LabelMode', 'name+percent+count');
```

`stats` contains `Category`, `Count`, and `Percent` and can be exported as a table.

## Web analysis interfaces

The web analyses are not identical to the standalone functions in this directory:

| Bridge | Input | Output |
| --- | --- | --- |
| `mdpEcFva` | ecModel identifier, target reaction, fraction of optimum | `fva_table`, `bar_chart` |
| `mdpKnockout` | ecModel identifier, a gene list or `all`, and a carbon-source reaction | `knockout_table`, `heatmap` |
| `mdpProteinAnalysis` | ecModel identifier and `subsystem`, `compartment`, or `gene` grouping | `usage_table`, `pie_data` |

Important methodological details are:

- `mdpEcFva` uses a default fraction of 0.9 and attempts a serial loop if parallel `ecFVA` fails.
- `mdpKnockout` reports `growth_ratio`, defined as the absolute post-deletion biomass flux divided by the wild-type biomass flux.
- `mdpProteinAnalysis` groups kcat weights as a proxy for enzyme demand. It is not a rigorous proteomics quantification based on solved enzyme-usage fluxes.

Publications and reports should identify the exact entry point and statistical definition used.

## Live Scripts

```matlab
setup
open('scripts/Analysis/chemostat_batch_simulation.mlx')
```

Live Scripts may contain project-specific model names, reaction identifiers, and file paths. Before execution, inspect the input cells for:

- the project root and model artifacts;
- biomass, carbon-source, and product reaction identifiers;
- experimental-data paths; and
- output directories for figures and tables.

## Output management

Store analysis results under:

```text
projects/<project>/analysis/
```

Record the model artifact, parameters, and reaction identifiers used for each result. For computationally intensive procedures such as FVA, retain both the numerical output and the plotting code rather than only the rendered figure.

## Troubleshooting

### `Undefined function 'progressbar'`

Confirm that RAVEN or the relevant utility directory is on the MATLAB path. Alternatively, use the web bridge, which attempts a serial implementation when parallel FVA fails.

### ecFVA mapping fails

Verify that the ecModel and GEM originate from the same network and that `mapRxnsToConv` recognizes their reaction identifiers. If the two models were renamed or modified independently, construct an explicit reaction mapping before comparison.

### The parallel pool cannot be started

The direct `ecFVA` implementation calls `parpool` and `parfor`. Confirm that Parallel Computing Toolbox is available, use the web serial fallback, or implement a dedicated serial loop for the local analysis.

### Empty groups, NaN values, or invalid logarithmic axes

Remove `NaN` and `Inf` values before plotting, and ensure that log10 ECDF inputs are positive. `plotEcFVA` ignores approximately zero-flux reactions, but multiple models must still share the same row ordering.
