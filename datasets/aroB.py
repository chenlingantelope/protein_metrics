# data is copied from the entropy project in data repo
# set working directory
import os
os.chdir("/")
# import seqIO to read fasta files
from Bio import SeqIO
from Bio.Seq import Seq
import os
from hmm_profile import reader
from esm.inverse_folding.util import load_coords


# read aroB_design.fasta using biopython
def load_dataset(dataset_name, recompute=False):
    if dataset_name =='aroB':
        # convert design_seqs and pdb sequence to the MSA coordinate
        # first create a gene profile using hmmbuild
        # check if the output file exists
        if (not os.path.exists("../seqs/aroB.hmm")) or recompute:
            # build a gene profile using hmmbuild
            assert os.system("hmmbuild -n aroB -o seqs/aroB.hmmbuild.out seqs/aroB.hmm seqs/aroB_msa.fasta") == 0
        with open("../seqs/aroB.hmm") as f:
            hmm = reader.read_single(f)
            # then align all three sets of fasta sequences to the gene profile so that they share the same coordinate system
        pdb = 'seqs/5eks.pdb'
        coords, pdb_seq = load_coords(pdb, "A")
        # write pdb_seq to a fasta file

        if (not os.path.exists("../seqs/aroB_design.a2m")) or recompute:
            assert os.system("hmmalign --outformat a2m -o seqs/aroB_design.a2m seqs/aroB.hmm seqs/aroB_design.fasta") == 0
        if (not os.path.exists("../seqs/5eks.a2m")) or recompute:
            with open("../seqs/5eks.fasta", "w") as f:
                f.write(f">5eks\n{pdb_seq}")
            f.close()
            assert os.system("hmmalign --outformat a2m -o seqs/5eks.a2m seqs/aroB.hmm seqs/5eks.fasta") == 0
        if (not os.path.exists("../seqs/aroB.a2m")) or recompute:
            assert os.system("hmmalign --outformat a2m -o seqs/aroB.a2m seqs/aroB.hmm seqs/aroB.fasta") == 0
        # read the aligned sequences
        design_seqs = list(SeqIO.parse("seqs/aroB_design.a2m", "fasta"))
        pdb_seq = list(SeqIO.parse("seqs/5eks.a2m", "fasta"))[0]
        wildtype_seq = list(SeqIO.parse("seqs/aroB.a2m", "fasta"))[0]
    else:
        raise ValueError("dataset_name should be either 'aroB' or 'infA'")
    # convert wildtype to SeqRecord object
    return design_seqs, pdb_seq, wildtype_seq, hmm, coords

