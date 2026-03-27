# ECOMAP Reconstruction GUI

Graphical user interface for enzyme-constraint metabolic model (ecGEM) reconstruction.

## Directory Structure

```
scripts/GUI/
├── Core/
│   ├── ReconstructionApp.m    # Main GUI application class
│   └── DockerChecker.m       # Docker installation/running checker
├── Utils/
│   ├── ModelStats.m          # Model statistics visualization
│   └── FilterReport.m        # Reaction filter report visualization
├── launchReconstructionGUI.m # Launcher function
└── README.md
```

## Quick Start

```matlab
% Add GUI to path and launch
cd scripts/GUI
launchReconstructionGUI()
```

## Features

### Workflow Steps

1. **Load Model** - Load metabolic model from SBML, JSON, or YAML
2. **Configure Parameters** - Set sigma, Ptot, f (protein constraint parameters)
3. **Convert to ecModel** - Build enzyme-constrained model
4. **Set Kcat Source** - Load existing kcat data or run prediction (requires Docker)
5. **Save Model** - Export reconstructed ecModel

### Parameter Descriptions

| Parameter | Description | Default |
|-----------|-------------|---------|
| sigma | Spheroid coverage factor | 0.5 |
| Ptot | Total protein content (g/gDW) | 0.3 |
| f | Fraction of active enzyme pool | 0.5 |

### ecModel Types

- **basic**: Irreversible split + protein pool bookkeeping
- **isozyme**: + Network isozyme expansion
- **integrated**: Enzyme pseudo-metabolites + usage reactions (recommended)

## Docker Integration

Kcat prediction tools (CatPred, DLKcat, UniKP) require Docker:

```matlab
% Check Docker status
result = DockerChecker.checkDocker();
disp(result);
```

Requirements:
- Docker Desktop installed
- Docker Desktop running (daemon active)
- Image: `hanhanxioyoge/kcat_prediction_catpred:v1.0` (for CatPred)

## Academic Style Guidelines

The GUI follows a minimalist academic design:
- Clean white/gray color scheme
- Clear step-by-step navigation
- Essential information only
- Professional typography (system fonts)

## Adding to MATLAB Path

For persistent path setup, add to `startup.m`:

```matlab
addpath(genpath('D:/project/ECOMAP/ECOMAP/scripts/GUI'));
```

## Troubleshooting

### GUI doesn't launch
- Ensure MATLAB R2022a or later
- Check COBRA Toolbox is installed

### Docker not found
- Install Docker Desktop
- Start Docker Desktop before prediction
- Click "Check Status" to verify connection

### Model conversion fails
- Verify model file format (SBML recommended)
- Check model contains gene-reaction associations
- Ensure UniProt data is available

## Code Mode

To switch to code-based workflow, use the underlying functions directly:

```matlab
% Load model
model = loadModel('iML1515.xml', 'TRADITION', modelDir, params);

% Convert to ecModel
ecModel = convertecModel(model, 'integrated', params);

% Save model
save('ecModel.mat', 'ecModel');
```
