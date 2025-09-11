# ECOMAP - ecModel Reconstruction Workflow 🧫

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![MATLAB R2022a+](https://img.shields.io/badge/MATLAB-R2022a%2B-blue.svg)](https://www.mathworks.com/)
[![COBRA Toolbox](https://img.shields.io/badge/COBRA_Toolbox-3.1%2B-orange.svg)](https://opencobra.github.io/cobratoolbox/)

> Construct ecModel from conventional GEM | Import ecModel developed by alternative methodologies

## 🌐 Environmental preparation
### Toolkit
- [COBRA Toolbox](https://github.com/opencobra/cobratoolbox) 3.1+
- [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) 2.9.2+
- [Gurobi Optimizer](https://www.gurobi.com/) (Or other MATLAB-supported solvers)

## Construct ecModel from conventional GEM

The construction of an enzyme-constrained metabolic model (ecModel) from a conventional genome-scale metabolic model (GEM) can be efficiently accomplished using the ECOMAP framework. This process involves the integration of enzyme capacity constraints, primarily through the incorporation of ***k*<sub>cat</sub> values** and **total protein constraints**, to improve the predictive accuracy of metabolic fluxes. Here we take `iML1515` as an example to build its ecModel:

1. Load the traditional COBRA model through the ``loadModel`` function

    ```
    Path = 'Select a path to store the model and related information';
    filename = 'iML1515.xml';
    fullFilePath = fullfile(Path, filename);
    model = loadModel(fullFilePath, 'COBRA');
    ```
2. Initialization parameters(``'ptot'`` represents the total protein content[g/gDW], ``'f'`` represents mass fraction of enzymes, ``'sigma'`` represents average enzyme saturation factor.)

    ```
    model = Parameterization(model, 'organism', 'Escherichia coli', ...
                                    'taxonomicid', 83333, ...
                                    'keggid', 'eco', ...
                                    'uniprottype', 'proteome', ...
                                    'uniprotid', 'UP000000625', ...
                                    'ptot', 0.55, ...
                                    'f', 0.55, ...
                                    'sigma', 0.5);
    ```
3. Converted to ecModel structure (Two types: `'complex'` and `'simple'`)

    ```
    ecModel = convertEcModel(model, 'complex');
    ```

4. Fill the ecModel with protein `uniprot_id`, `sequences`, and `relative molecular mass`

    ```
    ecModel = fillEnzymeInformation(ecModel);
    ```

5. Initially retrieve enzyme ***k*<sub>cat</sub>** through the `BRENDA` database

    ```
    ecModel = getkcatfromDatabase(ecModel);
    ```
    
6. Get metabolite SMILES format

    ```
    [ecModel, noSMILES] = findMetSmiles(ecModel)
    ```
 
7. Run DLKcat to predict vacant ***k*<sub>cat</sub>** values
   
    ```
    runDLKcat();
    ```

## Import ecModel developed by alternative methodologies

1. Load ecModel through the loadModel function
    ```
    Path = 'Select a path to store the model and related information';
    filename = 'ecYeastGEM.yml';
    fullFilePath = fullfile(Path, filename);
    model = loadModel(fullFilePath, 'GECKO');
    ```
2. Initialization parameters
    ```
    model = Parameterization(model, 'organism', 'Escherichia coli', ...
                                    'taxonomicid', 83333, ...
                                    'keggid', 'eco', ...
                                    'uniprottype', 'proteome', ...
                                    'uniprotid', 'UP000000625', ...
                                    'ptot', 0.55, ...
                                    'f', 0.55, ...
                                    'sigma', 0.5);
    ```
