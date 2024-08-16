import os
import pandas as pd
from data_loader import Dataset
from biotite.sequence import ProteinSequence, NucleotideSequence
from biotite.sequence.io import fasta
from biotite.application.blast import BlastWebApp
from biotite.database import entrez
from tempfile import gettempdir

# def mutate_wt(wt, x):
#     if x == '_wt':
#         return wt
#     else:
#         x = x.split('(')[1].split(')')[0]
#     if 'ins' in x:
#         pos = int(x.split('_')[0][1:])
#         ins = x.split('ins')[1]
#         assert wt[pos-1] == x[0]
#         return wt[:pos-1] + ProteinSequence(ins) + wt[pos-1:]
#     elif 'del' in x:
#         x = x.replace('del', '')
#         if '_' in x:
#             x1 = int(x.split('_')[0][1:])
#             x2 = int(x.split('_')[1][1:])
#             assert wt[x1-1] == x.split('_')[0][0]
#             assert wt[x2-1] == x.split('_')[1][0]
#             return wt[:x1-1] + wt[x2-1:]
#         else:
#             x1 = int(x[1:])
#             assert wt[x1-1] == x[0]
#             return wt[:x1-1] + wt[x1:]
#     else:
#         pos = int(x[1:-1])
#         assert wt[pos-1] == x[0]
#         return wt[:pos-1] + ProteinSequence(x[-1]) + wt[pos:]

dewachter2023 = Dataset('dewachter2023')
# urls = ['https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-35940-3/MediaObjects/41467_2023_35940_MOESM7_ESM.csv',
#         'https://raw.githubusercontent.com/odcambc/DIMPLE/master/tests/Kir.fa']

#Read the fabZ.fa folder to get the sequence
dewachter2023.get_wt('', 'fabZ', read_from_file=True, filename='fabZ.fa')

# url_files = ['score.txt', 'wt.fa']

# for i,x in enumerate(urls):
#     dewachter2023.download_from_url(x, url_files[i])
#
# # read the wt.fa file
# seq = fasta.FastaFile.read(os.path.join(dimple.dir, url_files[1]))
# seq = fasta.get_sequence(seq)
# wt = seq[582:1893].translate(complete=True)
# dimple.get_wt(wt, 'Kir21')
# # write the wt sequence to a file
# with open(os.path.join(dimple.dir, 'Kir21.fa'), 'w') as f:
#     f.write(f'>Kir21\n{wt}\n')
# f.close()
# # mutate the wildtype sequence
# data = pd.read_csv(os.path.join(dimple.dir, url_files[1]), sep='\t')
#
#
# variants = [mutate_wt(wt, x) for x in data['hgvs']]
# variantid = list(data['hgvs'])
# dimple.get_design_sequence(variants, variantid, data['score'])
#
# # get msa alignment
# dimple.get_msa(os.path.join(dimple.dir, 'wt.fa'))
# dimple.get_hmm()
#
# dimple.get_pdb_sequence('3SPI')
#
#
# from external.ProteinMPNN.protein_mpnn_utils import parse_PDB, ProteinMPNN
# import torch
# pdb_dict_list = parse_PDB(dimple.pdb_file)
# checkpoint = torch.load('external/ProteinMPNN/vanilla_model_weights/v_48_030.pt')
# noise_level_print = checkpoint['noise_level']
# hidden_dim = 128
# num_layers = 3
# model = ProteinMPNN(ca_only=False, num_letters=21, node_features=hidden_dim, edge_features=hidden_dim,
#                     hidden_dim=hidden_dim, num_encoder_layers=num_layers, num_decoder_layers=num_layers,
#                     k_neighbors=checkpoint['num_edges'])
# model.load_state_dict(checkpoint['model_state_dict'])
# model.eval()
#
# # align the pdb sequence to the wildtype sequence
#
# pdb_dict_list[0]['num_of_chains']
# pdb_dict_list[0]['seq']
# # biotite pairwise align pdb sequence to wildtype sequence
# from biotite.sequence.align import align_optimal, SubstitutionMatrix
#
# alignment = align_optimal(ProteinSequence(pdb_dict_list[0]['seq']), dimple.wt,
#                           matrix=SubstitutionMatrix.std_protein_matrix(), local=True)
#
# for a in alignment:
#     print(a)
#
# len(alignment[0])
# len(pdb_dict_list[0]['seq'])
# len(dimple.wt)
# alignment[0].trace
# # subset the pdb_dict coordinates to only residues that are aligned to the wildtype sequence
# temp = pdb_dict_list[0]['coords_chain_A']
# temp.keys()
# pdb_dict_list[0]['seq_chain_A']