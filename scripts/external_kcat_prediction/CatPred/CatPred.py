"""
Enzyme Kinetics Parameter Prediction Script

This script predicts enzyme kinetics parameters (kcat, Km, or Ki) using a pre-trained model.
It processes input data, generates predictions, and saves the results.

Usage:
    python demo_run.py --parameter <kcat|km|ki> --input_file <path_to_input_csv> --checkpoint_dir <path_to_pretrained_checkpoint_dir> [--use_gpu]

Dependencies:
    pandas, numpy, rdkit, IPython, argparse
"""

import time
import os, stat
import pandas as pd
import numpy as np
from IPython.display import Image, display
from rdkit import Chem
from IPython.display import display, Latex, Math
import argparse

ESM_MAX_LEN = 1022
VALID_AAS = set('ACDEFGHIKLMNPQRSTVWY')

def create_csv_sh(parameter, input_file_path, checkpoint_dir):
    df = pd.read_csv(input_file_path)

    df.columns = [c.strip() for c in df.columns]
    if not {'SMILES', 'sequence'}.issubset(df.columns):
        raise ValueError(f"Required columns are missing: { {'SMILES','sequence'} - set(df.columns) }")

    df['SMILES'] = df['SMILES'].astype(str).str.strip()
    df['sequence'] = df['sequence'].astype(str).str.strip().str.upper()

    keep_idx = []
    cleaned_smiles = {}

    for idx, row in df.iterrows():
        smi = row['SMILES']
        seq = row['sequence']

        # Sequence verification
        if not set(seq).issubset(VALID_AAS):
            continue
        if len(seq) == 0 or len(seq) > ESM_MAX_LEN:
            continue

        # SMILES parsing & normalizing
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            continue

        smi_can = Chem.MolToSmiles(mol)
        if parameter == 'kcat' and '.' in smi_can:
            smi_can = '.'.join(sorted(smi_can.split('.')))  # Fragment sorting to ensure certainty

        keep_idx.append(idx)
        cleaned_smiles[idx] = smi_can

    # Only legitimate lines are kept and normalized SMILES are written
    df_valid = df.loc[keep_idx].copy()
    df_valid['SMILES'] = [cleaned_smiles[i] for i in keep_idx]

    # Write cleaned input file (no index)
    base, _ = os.path.splitext(input_file_path)
    input_file_new_path = f"{base}_input.csv"
    df_valid.to_csv(input_file_new_path, index=False)

    with open('predict.sh', 'w') as f:
        f.write(f'''
        TEST_FILE_PREFIX={input_file_new_path[:-4]}
        RECORDS_FILE=${{TEST_FILE_PREFIX}}.json
        CHECKPOINT_DIR={checkpoint_dir}
        
        python ./scripts/create_pdbrecords.py --data_file ${{TEST_FILE_PREFIX}}.csv --out_file ${{RECORDS_FILE}}
        python predict.py --test_path ${{TEST_FILE_PREFIX}}.csv --preds_path ${{TEST_FILE_PREFIX}}_output.csv --checkpoint_dir $CHECKPOINT_DIR --uncertainty_method mve --smiles_column SMILES --individual_ensemble_predictions --protein_records_path $RECORDS_FILE
        ''')

    with open('predict.sh', 'rb') as _f:
        _data = _f.read().replace(b'\r\n', b'\n')
    with open('predict.sh', 'wb') as _f:
        _f.write(_data)

    st = os.stat('predict.sh')
    os.chmod('predict.sh', st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    
    return input_file_new_path[:-4]+'_output.csv'

def get_predictions(parameter, outfile):
    """
    Process prediction results and add additional metrics.

    Args:
        parameter (str): The kinetics parameter that was predicted.
        outfile (str): Path to the output CSV file from the prediction.

    Returns:
        pandas.DataFrame: Processed predictions with additional metrics.
    """
    df = pd.read_csv(outfile)
    pred_col, pred_logcol, pred_sd_totcol, pred_sd_aleacol, pred_sd_epicol = [], [], [], [], []

    unit = 'mM'
    if parameter == 'kcat':
        target_col = 'log10kcat_max'
        unit = 's^(-1)'
    elif parameter == 'km':
        target_col = 'log10km_mean'
    else:
        target_col = 'log10ki_mean'

    unc_col = f'{target_col}_mve_uncal_var'

    for _, row in df.iterrows():
        model_cols = [col for col in row.index if col.startswith(target_col) and 'model_' in col]

        unc = row[unc_col]
        prediction = row[target_col]
        prediction_linear = np.power(10, prediction)

        model_outs = np.array([row[col] for col in model_cols])
        epi_unc = np.var(model_outs)
        alea_unc = unc - epi_unc
        epi_unc = np.sqrt(epi_unc)
        alea_unc = np.sqrt(alea_unc)
        unc = np.sqrt(unc)

        pred_col.append(prediction_linear)
        pred_logcol.append(prediction)
        pred_sd_totcol.append(unc)
        pred_sd_aleacol.append(alea_unc)
        pred_sd_epicol.append(epi_unc)

    df[f'Prediction_({unit})'] = pred_col
    df['Prediction_log10'] = pred_logcol
    df['SD_total'] = pred_sd_totcol
    df['SD_aleatoric'] = pred_sd_aleacol
    df['SD_epistemic'] = pred_sd_epicol

    return df

def main(args):
    print(os.getcwd())

    outfile = create_csv_sh(args.parameter, args.input_file, args.checkpoint_dir)
    if outfile is None:
        return

    print('Predicting.. This will take a while..')

    if args.use_gpu:
        print(args.use_gpu)
        os.system("export PROTEIN_EMBED_USE_CPU=0;./predict.sh")
    else:
        os.system("export PROTEIN_EMBED_USE_CPU=1;./predict.sh")

    output_final = get_predictions(args.parameter, outfile)
    output_final.to_csv(outfile, index=False)  # 直接写回原路径
    print('Output saved to', outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict enzyme kinetics parameters.")
    parser.add_argument("--parameter", type=str, choices=["kcat", "km", "ki"], required=True,
                        help="Kinetics parameter to predict (kcat, km, or ki)")
    parser.add_argument("--input_file", type=str, required=True,
                        help="Path to the input CSV file")
    parser.add_argument("--use_gpu", action="store_true",
                        help="Use GPU for prediction (default is CPU)")
    parser.add_argument("--checkpoint_dir", type=str, required=True,
                        help="Path to the model checkpoint directory")

    args = parser.parse_args()
    args.parameter = args.parameter.lower()

    main(args)
