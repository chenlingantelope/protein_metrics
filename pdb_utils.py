from biotite.structure import filter_amino_acids
from biotite.sequence import ProteinSequence
from datasets.infA import infA
from datasets.dimple import dimple
from biotite.structure.io import load_structure
pdb = load_structure(infA.pdb_file)
pdb = load_structure(dimple.pdb_file)


def get_protein_seq(pdb):
    if len(pdb.shape) == 2:
        # multiple conformers submitted in the pdb file
        # only keep the same conformer
        pdb = pdb[0]
    pdb_aa = pdb[filter_amino_acids(pdb)]
    pos = pdb_aa.get_annotation('res_id')
    res = pdb_aa.get_annotation('res_name')
    seq = {}
    for i, p in enumerate(pos):
        seq[p] = ProteinSequence().convert_letter_3to1(res[i])

    for i in range(max(pos)):
        if i not in pos:
            seq[i] = 'X'
    aa_seq = ''
    for i in range(max(pos)):
        aa_seq += seq[i]
    return aa_seq
