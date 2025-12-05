#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-10-23

import os, sys, math, csv
import pickle
import numpy as np
from rdkit import Chem
from collections import defaultdict
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.metrics import mean_squared_error,r2_score


def load_pickle(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)

fingerprint_dict = load_pickle('input/fingerprint_dict.pickle')
atom_dict = load_pickle('input/atom_dict.pickle')
bond_dict = load_pickle('input/bond_dict.pickle')
edge_dict = load_pickle('input/edge_dict.pickle')
word_dict = load_pickle('input/sequence_dict.pickle')

def split_sequence(sequence, ngram):
    sequence = '-' + sequence + '='
    # print(sequence)
    # words = [word_dict[sequence[i:i+ngram]] for i in range(len(sequence)-ngram+1)]

    words = list()
    for i in range(len(sequence)-ngram+1) :
        try :
            words.append(word_dict[sequence[i:i+ngram]])
        except :
            word_dict[sequence[i:i+ngram]] = 0
            words.append(word_dict[sequence[i:i+ngram]])

    return np.array(words)
    # return word_dict

def create_atoms(mol):
    """Create a list of atom (e.g., hydrogen and oxygen) IDs
    considering the aromaticity."""
    # atom_dict = defaultdict(lambda: len(atom_dict))
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    # print(atoms)
    for a in mol.GetAromaticAtoms():
        i = a.GetIdx()
        atoms[i] = (atoms[i], 'aromatic')
    atoms = [atom_dict[a] for a in atoms]
    # atoms = list()
    # for a in atoms :
    #     try: 
    #         atoms.append(atom_dict[a])
    #     except :
    #         atom_dict[a] = 0
    #         atoms.append(atom_dict[a])

    return np.array(atoms)

def create_ijbonddict(mol):
    """Create a dictionary, which each key is a node ID
    and each value is the tuples of its neighboring node
    and bond (e.g., single and double) IDs."""
    # bond_dict = defaultdict(lambda: len(bond_dict))
    i_jbond_dict = defaultdict(lambda: [])
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bond = bond_dict[str(b.GetBondType())]
        i_jbond_dict[i].append((j, bond))
        i_jbond_dict[j].append((i, bond))
    return i_jbond_dict

def extract_fingerprints(atoms, i_jbond_dict, radius):
    """Extract the r-radius subgraphs (i.e., fingerprints)
    from a molecular graph using Weisfeiler-Lehman algorithm."""

    # fingerprint_dict = defaultdict(lambda: len(fingerprint_dict))
    # edge_dict = defaultdict(lambda: len(edge_dict))

    if (len(atoms) == 1) or (radius == 0):
        fingerprints = [fingerprint_dict[a] for a in atoms]

    else:
        nodes = atoms
        i_jedge_dict = i_jbond_dict

        for _ in range(radius):

            """Update each node ID considering its neighboring nodes and edges
            (i.e., r-radius subgraphs or fingerprints)."""
            fingerprints = []
            for i, j_edge in i_jedge_dict.items():
                neighbors = [(nodes[j], edge) for j, edge in j_edge]
                fingerprint = (nodes[i], tuple(sorted(neighbors)))
                # fingerprints.append(fingerprint_dict[fingerprint])
                # fingerprints.append(fingerprint_dict.get(fingerprint))
                try :
                    fingerprints.append(fingerprint_dict[fingerprint])
                except :
                    fingerprint_dict[fingerprint] = 0
                    fingerprints.append(fingerprint_dict[fingerprint])

            nodes = fingerprints

            """Also update each edge ID considering two nodes
            on its both sides."""
            _i_jedge_dict = defaultdict(lambda: [])
            for i, j_edge in i_jedge_dict.items():
                for j, edge in j_edge:
                    both_side = tuple(sorted((nodes[i], nodes[j])))
                    # edge = edge_dict[(both_side, edge)]
                    # edge = edge_dict.get((both_side, edge))
                    try :
                        edge = edge_dict[(both_side, edge)]
                    except :
                        edge_dict[(both_side, edge)] = 0
                        edge = edge_dict[(both_side, edge)]

                    _i_jedge_dict[i].append((j, edge))
            i_jedge_dict = _i_jedge_dict

    return np.array(fingerprints)

def create_adjacency(mol):
    adjacency = Chem.GetAdjacencyMatrix(mol)
    return np.array(adjacency)

def dump_dictionary(dictionary, filename):
    with open(filename, 'wb') as file:
        pickle.dump(dict(dictionary), file)

def load_tensor(file_name, dtype):
    return [dtype(d).to(device) for d in np.load(file_name + '.npy', allow_pickle=True)]

class Predictor(object):
    def __init__(self, model):
        self.model = model

    def predict(self, data):
        predicted_value = self.model.forward(data)

        return predicted_value

class KcatPrediction(nn.Module):
    def __init__(self, device, n_fingerprint, n_word, dim, layer_gnn, window, layer_cnn, layer_output):
        super(KcatPrediction, self).__init__()
        self.embed_fingerprint = nn.Embedding(n_fingerprint, dim)
        self.embed_word = nn.Embedding(n_word, dim)
        self.W_gnn = nn.ModuleList([nn.Linear(dim, dim)
                                    for _ in range(layer_gnn)])
        self.W_cnn = nn.ModuleList([nn.Conv2d(
                     in_channels=1, out_channels=1, kernel_size=2*window+1,
                     stride=1, padding=window) for _ in range(layer_cnn)])
        self.W_attention = nn.Linear(dim, dim)
        self.W_out = nn.ModuleList([nn.Linear(2*dim, 2*dim)
                                    for _ in range(layer_output)])
        # self.W_interaction = nn.Linear(2*dim, 2)
        self.W_interaction = nn.Linear(2*dim, 1)

        self.device = device
        self.dim = dim
        self.layer_gnn = layer_gnn
        self.window = window
        self.layer_cnn = layer_cnn
        self.layer_output = layer_output

    def gnn(self, xs, A, layer):
        for i in range(layer):
            hs = torch.relu(self.W_gnn[i](xs))
            xs = xs + torch.matmul(A, hs)
        # return torch.unsqueeze(torch.sum(xs, 0), 0)
        return torch.unsqueeze(torch.mean(xs, 0), 0)

    def attention_cnn(self, x, xs, layer):
        """The attention mechanism is applied to the last layer of CNN."""

        xs = torch.unsqueeze(torch.unsqueeze(xs, 0), 0)
        for i in range(layer):
            xs = torch.relu(self.W_cnn[i](xs))
        xs = torch.squeeze(torch.squeeze(xs, 0), 0)

        h = torch.relu(self.W_attention(x))
        hs = torch.relu(self.W_attention(xs))
        weights = torch.tanh(F.linear(h, hs))
        ys = torch.t(weights) * hs

        # return torch.unsqueeze(torch.sum(ys, 0), 0)
        return torch.unsqueeze(torch.mean(ys, 0), 0)

    def forward(self, inputs):

        fingerprints, adjacency, words = inputs

        layer_gnn = 3
        layer_cnn = 3
        layer_output = 3

        """Compound vector with GNN."""
        fingerprint_vectors = self.embed_fingerprint(fingerprints)
        compound_vector = self.gnn(fingerprint_vectors, adjacency, layer_gnn)

        """Protein vector with attention-CNN."""
        word_vectors = self.embed_word(words)
        protein_vector = self.attention_cnn(compound_vector,
                                            word_vectors, layer_cnn)

        """Concatenate the above two vectors and output the interaction."""
        cat_vector = torch.cat((compound_vector, protein_vector), 1)
        for j in range(layer_output):
            cat_vector = torch.relu(self.W_out[j](cat_vector))
        interaction = self.W_interaction(cat_vector)
        # print(interaction)

        return interaction

    def __call__(self, data, train=True):

        inputs, correct_interaction = data[:-1], data[-1]
        predicted_interaction = self.forward(inputs)
        print(predicted_interaction)

        if train:
            loss = F.mse_loss(predicted_interaction, correct_interaction)
            return loss
        else:
            correct_values = correct_interaction.to('cpu').data.numpy()
            predicted_values = predicted_interaction.to('cpu').data.numpy()[0]
            # correct_values = np.concatenate(correct_values)
            # predicted_values = np.concatenate(predicted_values)
            # ys = F.softmax(predicted_interaction, 1).to('cpu').data.numpy()
            # predicted_values = list(map(lambda x: np.argmax(x), ys))
            print(correct_values)
            print(predicted_values)
            # predicted_scores = list(map(lambda x: x[1], ys))
            return correct_values, predicted_values

def _get_first(row, *keys):
    for k in keys:
        if k in row and row[k] is not None:
            v = str(row[k]).strip()
            if v != '':
                return v
    return ''

def main() :
    if len(sys.argv) < 3:
        print("Usage: python Dlkcat.py <input.csv> <output.csv>")
        sys.exit(1)

    inputfile = sys.argv[1:][0]
    outputfile = sys.argv[1:][1]
    #print(inputfile)

    if not os.access(inputfile, os.R_OK):
        print('DLKcat cannot find the input file ' + inputfile)
        sys.exit(1)

    fingerprint_dict = load_pickle('input/fingerprint_dict.pickle')
    atom_dict = load_pickle('input/atom_dict.pickle')
    bond_dict = load_pickle('input/bond_dict.pickle')
    word_dict = load_pickle('input/sequence_dict.pickle')
    n_fingerprint = len(fingerprint_dict)
    n_word = len(word_dict)

    radius=2
    ngram=3

    dim=10
    layer_gnn=3
    side=5
    window=11
    layer_cnn=3
    layer_output=3
    lr=1e-3
    lr_decay=0.5
    decay_interval=10
    weight_decay=1e-6
    iteration=100

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # torch.manual_seed(1234)
    Kcat_model = KcatPrediction(device, n_fingerprint, n_word, 2*dim, layer_gnn, window, layer_cnn, layer_output).to(device)
    Kcat_model.load_state_dict(torch.load('input/all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50', map_location=device))
    # print(state_dict.keys())
    # model.eval()
    predictor = Predictor(Kcat_model)

    out_header = [
        'ReactionName', 'Organism', 'GeneID', 'ProteinID', 'EC Number', 'MetaNetXID',
        'Substrate_name', 'Substrate_smiles', 'InChIKey','Protein_sequence',
        'predicted_kcat'
    ]
    
    with open(inputfile, 'r', encoding='utf-8-sig', newline='') as infh, \
        open(outputfile, 'w', encoding='utf-8', newline='') as outfh:

        sample = infh.read(4096); infh.seek(0)

        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=',\t;')
        except csv.Error:
            dialect = csv.excel

        reader = csv.DictReader(infh, dialect=dialect)
        writer = csv.DictWriter(outfh, fieldnames=out_header)
        writer.writeheader()    

        for row in reader:
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

            pred_kcat = 'None'
            try:
                if smiles and "." not in smiles:
                    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                    if mol is not None:
                        atoms      = create_atoms(mol)
                        ijbond     = create_ijbonddict(mol)
                        fps        = extract_fingerprints(atoms, ijbond, radius)
                        adj        = create_adjacency(mol)
                        words      = split_sequence(seq, ngram)

                        fps_t  = torch.LongTensor(fps).to(device)
                        adj_t  = torch.FloatTensor(adj).to(device)
                        words_t= torch.LongTensor(words).to(device)

                        inputs  = [fps_t, adj_t, words_t]
                        pred    = predictor.predict(inputs)
                        log2val = pred.item()
                        pred_kcat = f"{math.pow(2, log2val):.4f}"
            except Exception:
                pred_kcat = 'None'

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
                'predicted_kcat'   : pred_kcat
            })

if __name__ == '__main__' :
    main()
