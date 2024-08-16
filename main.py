import pandas as pd
import torch

from datasets.aroB import load_dataset
from basic_metrics import all_basic_metrics
from msa_metrics import hmm_likelihood

from datasets.infA import infA
from datasets.aroB import aroB
from datasets.dimple import dimple
infA.design_seqs_score
aroB.pdb_file
design_seqs, pdb_seq, wildtype_seq, hmm, coords = load_dataset('aroB')
# Check that the sequence length correspond to the gene profile length
for seq in design_seqs + [pdb_seq] + [wildtype_seq]:
    n_upper = sum([i.isupper() for i in seq.seq])
    n_gaps = sum([i == '-' for i in seq.seq])
    assert n_upper + n_gaps == len(hmm.steps)

basic_metrics = all_basic_metrics(design_seqs, wildtype_seq)

# HMM likelihood of sequence to HMM alignment
hmm_likelihood = pd.Series([hmm_likelihood(seq.seq, hmm) for seq in design_seqs],
                           index=[seq.id for seq in design_seqs])


from structure_metrics import esm_if_score
esm_ll = esm_if_score(pdb_seq, design_seqs[0], coords)


from external.ProteinMPNN.protein_mpnn_utils import parse_PDB, ProteinMPNN

pdb_dict_list = parse_PDB('seqs/5eks.pdb')
checkpoint = torch.load('external/ProteinMPNN/vanilla_model_weights/v_48_030.pt')
noise_level_print = checkpoint['noise_level']
hidden_dim = 128
num_layers = 3
model = ProteinMPNN(ca_only=False, num_letters=21, node_features=hidden_dim, edge_features=hidden_dim,
                    hidden_dim=hidden_dim, num_encoder_layers=num_layers, num_decoder_layers=num_layers,
                    k_neighbors=checkpoint['num_edges'])
model.load_state_dict(checkpoint['model_state_dict'])
model.eval()