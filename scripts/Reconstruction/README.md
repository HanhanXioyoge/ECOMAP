# ECOMAP Reconstruction Module

This module normalizes conventional GEMs or external ecGEMs into the ECOMAP model representation and supports enzyme-constraint topology construction, protein and metabolite annotation, kcat prediction and matching, and model export. The numerical routines invoked by the web reconstruction workflow are also located in this directory.

Return to the [project overview](../../README.md).

## Supported inputs and outputs

### Input model types

`loadModel(filename, modeltype, modelDir, parameters)` accepts the following `modeltype` values:

- `TRADITION`: conventional GEM
- `ECOMAP`: ECOMAP ecModel
- `SMOMENT`: sMOMENT model
- `ECMPY`: ECMpy model
- `GECKO`: GECKO model

Model-type matching is case-insensitive. The current loader supports `.xml`, `.json`, `.yml`, and `.yaml`; it does not read `.mat` files.

### ecModel topologies

`convertecModel(model, ecModeltype, parameters)` supports:

| Topology | Representation |
| --- | --- |
| `basic` | Irreversible reaction splitting and total-protein-pool bookkeeping |
| `isozyme` | The basic representation with network-level isozyme expansion |
| `integrated` | Enzyme pseudo-metabolites, enzyme-usage reactions, and the integrated coupling structure |

For conventional GEM conversion, `convertecModel` expects `model.type == 'TRADITION'`.

## Dependencies

- Run `setup` from the repository root before using the module.
- RAVEN Toolbox, COBRA Toolbox, and a functional solver must be available on the MATLAB path.
- UniProt, Complex Portal, PubChem, ChEBI, and MetaNetX retrieval requires network access unless a local cache is available.
- CatPred, DLKcat, and UniKP have independent Docker, model-asset, and Python requirements under `scripts/external_kcat_prediction/`.

## Recommended workflow

### 1. Load project parameters

```matlab
root = fileparts(which('setup'));
manager = fullfile(root, 'projects', 'eciML1515', ...
    'eciML1515ParameterManagement.m');
params = ParameterManager.getParams(manager);

% Reconstruct host-specific paths after moving or cloning the project.
params.projectDir        = fullfile(root, 'projects', 'eciML1515');
params.path              = params.projectDir;
params.projectJson       = fullfile(params.projectDir, 'project.json');
params.parameterManager  = manager;
params.modelsDir         = fullfile(params.projectDir, 'models');
params.reconstructionDir = fullfile(params.projectDir, 'reconstruction');
params.analysisDir       = fullfile(params.projectDir, 'analysis');
```

At minimum, verify the following parameters:

- `InitialModel` and `modeltype`
- `sigma`, `Ptot`, and `f`
- `org_name` and `taxonomicID`
- `uniprot.type`, `uniprot.ID`, `uniprot.geneIDfield`, and `uniprot.reviewed`
- `c_source` and `bioRxn`

### 2. Import and normalize the GEM

```matlab
gem = loadModel(params.InitialModel, params.modeltype, ...
    params.modelsDir, params);
```

The loader creates the common model fields, MIRIAM containers, and a model identifier. A relative filename is resolved against `modelsDir`.

### 3. Construct ecModel topologies

```matlab
ecBasic = convertecModel(gem, 'basic', params);
ecIso   = convertecModel(gem, 'isozyme', params);
ecInt   = convertecModel(gem, 'integrated', params);
```

During conversion, the implementation identifies pure intercompartmental transport reactions, excludes them from enzyme constraints, reads or downloads `reconstruction/uniprot.tsv`, and constructs the `enzymeConstraints` structure.

### 4. Add annotations

Use the following entry points as required by the project:

| Function | Purpose | Principal local artifact |
| --- | --- | --- |
| `getECnumber` | Map EC numbers from UniProt annotations | `reconstruction/uniprot.tsv` |
| `getComplexdata` | Retrieve and cache Complex Portal records | `reconstruction/ComplexPortal.json` |
| `applyComplexdata` | Apply curated or proposed complex composition to the model | Updated model structure |
| `addMetMetaNetXID` | Add MetaNetX metabolite cross-references | Updated metabolite annotations |
| `getMetinfo` | Add SMILES and InChIKey values from local and remote sources | `reconstruction/metInfo.tsv` |

`applyComplexdata` requires `uniprot.tsv` in `parameters.reconstructionDir`. Remote services may impose rate limits or return incomplete records; retain and examine the unmatched outputs returned by the annotation routines.

### 5. Predict and reconcile kcat values

Principal entry points are:

| Function | Purpose |
| --- | --- |
| `writeInputFile` | Generate sequence and substrate inputs for external predictors |
| `getPrediction` | Read results from a selected deep-learning predictor |
| `dbKcatMatch` | Align predictions with BRENDA/SABIO data |
| `fuzzyKcatMatch` | Match by EC number, substrate, organism, and phylogenetic distance |
| `completeKcatMatch` | Aggregate matches from multiple predictors and databases |
| `fillCustomKcats` | Apply manually curated values from `kcatData/customKcats.csv` |
| `fillMissingKcatWithMedian` | Impute missing values with median estimates |
| `mergeKcats` | Reconcile sources before values are written to the model |

Prediction CSV files are normally stored under `reconstruction/kcatData/`. `getPrediction` reads generated results; it does not install or train the external prediction models.

### 6. Save model artifacts

Use `scripts/utilities/saveModel.m` to export model structures. Store the three topologies in the project `models/` directory and use names recognized by the web project-state scanner:

```text
<project>-basic.mat
<project>-isozyme.mat
<project>-integrated.mat
```

## Web reconstruction sequence

The web project page infers reconstruction status from the following sequence and its artifacts:

```text
initialize → loadGem → convert → annotate → predictKcat
           → compareKcat → mergeKcat → growthSave
```

The corresponding bridges are:

| Stage | Bridge |
| --- | --- |
| Initialization and import | `mdpInitProject`, `mdpLoadModel` |
| Conversion | `mdpConvertecModel` |
| Annotation | `mdpAnnotate` |
| Prediction | `mdpDlPredict` |
| Comparison and merging | `mdpKcatCompare`, `mdpKcatMerge` |
| Growth validation | `mdpGrowthPredict` |

Bridge model identifiers are registered in the current MATLAB Engine process. Models must be reloaded after the backend is restarted.

## Principal files

```text
Reconstruction/
├── loadModel.m                    # Model import and normalization
├── convertecModel.m               # GEM-to-ecModel topology construction
├── buildEnzConstrRxnSet.m         # Enzyme-constrained reaction selection
├── getECnumber.m                  # EC-number annotation
├── getMetinfo.m                   # Metabolite annotation
├── getComplexdata.m               # Complex-data retrieval and caching
├── applyComplexdata.m             # Complex composition integration
├── writeInputFile.m               # Predictor input generation
├── getPrediction.m                # Prediction-result import
├── dbKcatMatch.m                  # Database matching
├── completeKcatMatch.m            # Multi-source match aggregation
├── fillCustomKcats.m              # Curated kcat integration
└── mergeKcats.m                   # kcat reconciliation
```

## Troubleshooting

### `ParameterManager is not set`

Pass the absolute project parameter-manager path to `ParameterManager.getParams`, or pass a complete `params` structure explicitly to each function.

### `Unsupported file extension: .mat`

This is the expected behavior of the reconstruction loader. Use an XML, JSON, or YAML source model. If the MATLAB variable structure is already known, load a `.mat` file with MATLAB `load` and use only downstream routines compatible with that structure.

### `uniprot.tsv` or `ComplexPortal.json` cannot be found

Confirm that `params.reconstructionDir` identifies the active project `reconstruction/` directory. Initial retrieval requires network access; retained files can serve as offline caches.

### Conversion produces no enzyme-constrained reactions

Verify that the GEM contains valid `genes`, `grRules` or `rules`, reaction–gene associations, and EC/UniProt annotations. Pure transport reactions are deliberately excluded.

### kcat results are not merged

Check predictor names, CSV columns, reaction and protein identifiers, and the `reconstruction/kcatData/` path. Complete the match analysis before calling `completeKcatMatch` and `mergeKcats`; curated values should be placed in `customKcats.csv`.
