from copy import deepcopy
import numpy as np
from biotite.sequence.io.fasta import FastaFile

def msa_to_unmapped_seq_coord(seq, reverse=False):
    """
    Create a mapping from the MSA coordinate to the original sequence coordinate or vice versa
    :param seq: str, sequence in A2M format
    :param reverse: bool, if True, create a mapping from the original sequence coordinate to the MSA coordinate
    """

    j = -1
    k = -1
    mapping = {}
    for i in range(len(seq)):
        if seq[i].isupper():
            j+=1
            k+=1
        elif seq[i] == '-':
            j+=1
        elif seq[i].islower():
            k+=1
        if reverse:
            mapping[k] = j
        else:
            mapping[j] = k
    return mapping


def map_coords(coords, design_seq, pdb_seq):
    """
    Map the coordinates from the MSA to the protein sequence with structural information
    :param coords: list of np.array, coordinates of the protein structure
    :param design_seq: sequence to evaluate
    :param pdb_seq: sequence with structural information
    :return: coords, list of np.array, coordinates of the protein structure with the same length as the design_seq
    """
    mapped_coord = []
    masked_pos = deepcopy(coords[0])
    # replace number with np.inf
    masked_pos[:,:] = np.inf
    design_to_msa = msa_to_unmapped_seq_coord(design_seq, reverse=True)
    msa_to_pdb = msa_to_unmapped_seq_coord(pdb_seq)

    for i in range(len(design_to_msa)):
        msa_pos = design_to_msa[i]
        if msa_pos != -1:
            pdb_pos = msa_to_pdb[msa_pos]
            if pdb_pos != -1:
                mapped_coord.append(coords[pdb_pos])
            else:
                mapped_coord.append(masked_pos)
    return mapped_coord


def aligned_to_raw_seq(aligned_seq):
    """
    Convert aligned sequence to raw sequence
    :param aligned_seq: A2M format sequence
    :return: return the input sequence before aligning to HMM profile
    """
    return aligned_seq.replace("-", "").upper()


def read_dna(file, single=True):
    seqs = FastaFile.read(file)
    genes = {}
    for genename, seq in seqs.items():
        genes[genename] = seq
    if single:
        return seq
    else:
        return genes

