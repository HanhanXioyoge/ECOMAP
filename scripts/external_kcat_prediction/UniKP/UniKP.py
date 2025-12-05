#!/usr/bin/env python
# coding: utf-8
"""
UniKP.py — CSV → CSV inference, preserving original UniKP feature logic.

Usage:
    python UniKP.py <input.csv> <output.csv>

Input CSV must contain headers:
    ReactionName,Organism,GeneID,ProteinID,EC Number,Substrate_name,Substrate_smiles,Protein_sequence

Output CSV = input columns + one extra column:
    predicted_kcat

Notes:
- Keeps the original UniKP feature pipeline:
    * SMILES -> SMILES-Transformer (trfm_12_23000.pkl, vocab.pkl)  → 256-d embedding
    * Protein sequence -> ProtT5 (mean over tokens, drop last special) → H-d embedding
    * Concatenate [SMILES, Protein] -> regress with a pretrained sklearn regressor
- The pretrained regressor is loaded from 'unikp_regressor.pkl' (relative to script directory).
- By default we assume the regressor output is log10(kcat). Set ASSUME_LOG10=0 to bypass 10**y.
"""

import os
import sys
import csv
import math
import re
import gc
import pickle
from typing import List, Tuple

import numpy as np
import torch
from transformers import T5EncoderModel, T5Tokenizer

# SMILES-Transformer dependencies (from original UniKP codebase)
from build_vocab import WordVocab
from pretrain_trfm import TrfmSeq2seq
from utils import split as smiles_tokenize

# ----------------------- Configuration -----------------------
# If your regressor outputs linear kcat (not log10), set to 0.
ASSUME_LOG10 = int(os.environ.get("ASSUME_LOG10", "1"))

# Default resource files (relative to script dir unless absolute paths)
DEFAULT_VOCAB_PKL   = "vocab.pkl"
DEFAULT_TRFM_CKPT   = "trfm_12_23000.pkl"
DEFAULT_REGRESSORPKL= "UniKP_for_kcat.pkl"


# ----------------------- Path helpers ------------------------
def script_dir() -> str:
    return os.path.dirname(os.path.abspath(__file__))

def resolve_path(path_like: str) -> str:
    """Resolve relative paths against the script directory."""
    if os.path.isabs(path_like):
        return path_like
    return os.path.join(script_dir(), path_like)


# ----------------------- CSV helpers -------------------------
def _get_first(row, *keys) -> str:
    """Pick the first non-empty value from candidate keys in a CSV row."""
    for k in keys:
        if k in row and row[k] is not None:
            v = str(row[k]).strip()
            if v != "":
                return v
    return ""


# ----------------------- UniKP: SMILES encoder ----------------
def smiles_to_vec(Smiles: List[str]) -> np.ndarray:
    """
    Original UniKP behavior:
    - Tokenize with utils.split
    - Build padded ids with SOS/EOS, seq_len=220
    - Encode with TrfmSeq2seq.encode (transpose xid to time-major)
    """
    pad_index = 0
    unk_index = 1
    eos_index = 2
    sos_index = 3
    seq_len   = 220

    vocab = WordVocab.load_vocab(resolve_path(DEFAULT_VOCAB_PKL))

    def get_inputs(sm):
        sm = sm.split()
        if len(sm) > 218:
            # Keep the head+tail to cap at 220 incl. SOS/EOS (109 + 109)
            sm = sm[:109] + sm[-109:]
        ids = [vocab.stoi.get(token, unk_index) for token in sm]
        ids = [sos_index] + ids + [eos_index]
        seg = [1] * len(ids)
        padding = [pad_index] * (seq_len - len(ids))
        ids.extend(padding)
        seg.extend(padding)
        return ids, seg

    def get_array(smiles_tokens):
        x_id, x_seg = [], []
        for sm in smiles_tokens:
            a, b = get_inputs(sm)
            x_id.append(a)
            x_seg.append(b)
        return torch.tensor(x_id), torch.tensor(x_seg)

    # Build model and load weights (as in original code)
    trfm = TrfmSeq2seq(len(vocab), 256, len(vocab), 4)
    trfm.load_state_dict(torch.load(resolve_path(DEFAULT_TRFM_CKPT), map_location="cpu"))
    trfm.eval()

    # Tokenize smiles (string -> space-delimited tokens per original)
    x_split = [smiles_tokenize(sm) for sm in Smiles]
    xid, xseg = get_array(x_split)   # xid shape: [N, 220]
    X = trfm.encode(torch.t(xid))    # encode expects [seq_len, batch]; returns [batch, 256]
    return X                  


# ----------------------- UniKP: ProtT5 encoder ----------------
def Seq_to_vec(Sequence: List[str]) -> np.ndarray:
    """
    Original UniKP behavior with minor safety:
    - Truncate >1000 aa to 500 head + 500 tail
    - Replace U/Z/O/B with X (after spacing)
    - Use ProtT5 encoder; mean-pool over valid tokens (exclude last special)
    - Returns [N, hidden] array (float)
    """
    # Truncate very long sequences the same way as original code
    Sequence = list(Sequence)
    for i in range(len(Sequence)):
        if len(Sequence[i]) > 1000:
            Sequence[i] = Sequence[i][:500] + Sequence[i][-500:]

    # Build spaced sequences (e.g., "A B C ...")
    sequences_spaced = []
    for seq in Sequence:
        seq = seq or ""
        if not seq:
            sequences_spaced.append("")
            continue
        # original code builds a spaced string by iterating chars
        seq_spaced = " ".join(list(seq))
        seq_spaced = re.sub(r"[UZOB]", "X", seq_spaced)
        sequences_spaced.append(seq_spaced)

    CACHE     = os.path.abspath(__file__)
    tokenizer = T5Tokenizer.from_pretrained("prot_t5_xl_uniref50", cache_dir=CACHE, do_lower_case=False)
    model     = T5EncoderModel.from_pretrained("prot_t5_xl_uniref50", cache_dir=CACHE)

    print(torch.cuda.is_available())
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    model = model.to(device).eval()

    features = []
    for i in range(len(sequences_spaced)):
        print('For sequence ', str(i + 1))
        seq_i = sequences_spaced[i]
        ids = tokenizer.batch_encode_plus([seq_i], add_special_tokens=True, padding=True)
        input_ids     = torch.tensor(ids['input_ids']).to(device)
        attention_mask= torch.tensor(ids['attention_mask']).to(device)
        with torch.no_grad():
            embedding = model(input_ids=input_ids, attention_mask=attention_mask)
        embedding = embedding.last_hidden_state.cpu().numpy()  # [1, L, H]
        att = attention_mask.cpu().numpy()

        # Mean over valid tokens, drop last special token (match original)
        for seq_num in range(len(embedding)):
            seq_len = int((att[seq_num] == 1).sum())
            if seq_len <= 1:
                features.append(np.zeros(embedding.shape[-1], dtype=float))
            else:
                seq_emd = embedding[seq_num][:seq_len - 1]  # [L-1, H]
                features.append(seq_emd.mean(axis=0))

        gc.collect()

    features = np.asarray(features, dtype=float)  # [N, H]
    return features


# ----------------------- Regressor ----------------------------
def load_regressor(pkl_path: str):
    pkl = resolve_path(pkl_path)
    if not os.path.isfile(pkl):
        raise FileNotFoundError(
            f"Pretrained regressor not found: {pkl}\n"
            f"Please provide a sklearn regressor pickle named '{DEFAULT_REGRESSORPKL}' "
            f"in the script directory, or change DEFAULT_REGRESSORPKL."
        )
    with open(pkl, "rb") as f:
        reg = pickle.load(f)
    if not hasattr(reg, "predict"):
        raise TypeError("Loaded model has no .predict method.")
    return reg

def main():
    if len(sys.argv) < 3:
        print("Usage: python UniKP.py <input.csv> <output.csv>")
        sys.exit(1)

    inputfile  = sys.argv[1]
    outputfile = sys.argv[2]

    if not os.access(inputfile, os.R_OK):
        print('UniKP cannot find the input file ' + inputfile)
        sys.exit(1)

    # Read CSV with delimiter sniffing (compatible with ',', '\t', ';')
    with open(inputfile, 'r', encoding='utf-8-sig', newline='') as infh:
        sample = infh.read(4096); infh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=',\t;')
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(infh, dialect=dialect)
        rows = [r for r in reader]

    out_header = [
        'ReactionName', 'Organism', 'GeneID', 'ProteinID', 'EC Number', 'MetaNetXID',
        'Substrate_name', 'Substrate_smiles', 'InChIKey', 'Protein_sequence',
        'predicted_kcat'
    ]

    # Extract columns (robust to minor header variants)
    items = []
    for row in rows:
        rxn      = _get_first(row, 'ReactionName', 'Reaction', 'rxn', 'reaction_name')
        org      = _get_first(row, 'Organism', 'organism')
        gene     = _get_first(row, 'GeneID', 'Gene', 'gene_id', 'gene')
        protein  = _get_first(row, 'ProteinID', 'protein_id', 'Protein', 'protein')
        ecnum    = _get_first(row, 'EC Number', 'EC_Number', 'ECNumber', 'EC', 'ec_number')
        eataid   = _get_first(row, 'MetaNetXID')
        sub_name = _get_first(row, 'Substrate_name', 'Substrate', 'substrate_name', 'substrate')
        smiles   = _get_first(row, 'Substrate_smiles', 'SMILES', 'smiles')
        inchi_key= _get_first(row, 'InChIKey', 'InChI_Key')
        seq      = _get_first(row, 'Protein_sequence', 'Sequence', 'protein_sequence', 'sequence')
        items.append((rxn, org, gene, protein, ecnum, eataid, sub_name, smiles, inchi_key, seq))

    # Decide which rows to predict on: require non-empty SMILES without '.'.
    valid_idx = [i for i, it in enumerate(items) if it[7] and ('.' not in it[7])]
    print(f"[UniKP] Total rows: {len(items)}, valid for prediction: {len(valid_idx)}")

    pred_map = {}  # row_index -> predicted kcat (float)

    if len(valid_idx) > 0:
        # Build feature inputs per original UniKP logic
        smiles_list = [items[i][7] for i in valid_idx]
        seq_list    = [items[i][9] for i in valid_idx]

        # UniKP SMILES feature
        print(f"[UniKP] Encoding SMILES with SMILES-Transformer ...")
        # The original smiles_to_vec expects a list of "space tokenized" strings,
        # so we join tokens produced by utils.split with spaces.
        smiles_token_strs = [" ".join(smiles_tokenize(sm)) for sm in smiles_list]
        X_smiles = smiles_to_vec(smiles_token_strs)  # [n, 256]

        print(f"[UniKP] Encoding sequences with ProtT5 ...")
        X_seq = Seq_to_vec(seq_list)  # [n, H]

        X = np.concatenate([X_smiles, X_seq], axis=1)
        print(f"[UniKP] Feature shape = {X.shape}")

        reg = load_regressor(DEFAULT_REGRESSORPKL)
        yhat = reg.predict(X).astype(float)

        if ASSUME_LOG10:
            yhat = np.power(10.0, yhat)

        pred_map = {idx: float(v) for idx, v in zip(valid_idx, yhat)}

    with open(outputfile, 'w', encoding='utf-8', newline='') as outfh:
        writer = csv.DictWriter(outfh, fieldnames=out_header)
        writer.writeheader()
        for i, (rxn, org, gene, protein, ecnum, eataid, sub_name, smiles, inchi_key, seq) in enumerate(items):
            pk = pred_map.get(i, None)
            writer.writerow({
                'ReactionName'     : rxn,
                'Organism'         : org,
                'GeneID'           : gene,
                'ProteinID'        : protein,
                'EC Number'        : ecnum,
                'MetaNetXID'       : eataid,
                'Substrate_name'   : sub_name,     
                'Substrate_smiles' : smiles,       
                'InChIKey'         : inchi_key,
                'Protein_sequence' : seq,         
                'predicted_kcat'   : (f"{pk:.6f}" if isinstance(pk, float) else 'None')
            })

    print(f"[UniKP] Finished. Wrote predictions to: {outputfile}")


if __name__ == '__main__':
    main()
