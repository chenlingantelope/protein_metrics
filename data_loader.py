from biotite.sequence import ProteinSequence
from biotite.sequence.align import align_optimal
import subprocess
import requests
import os
from hmm import HMM
from constants import UNIREF90
from utils import read_dna

# define dataset class
class Dataset:
    def __init__(self, name, local=False):
        self.name = name
        self.dir = f"datasets/{name}"
        if not os.path.exists(self.dir):
            os.makedirs(self.dir)
        self.url, self.url_filename = [], []

    def __str__(self):
        return f"Dataset: {self.name}"

    def __repr__(self):
        return f"Dataset: {self.name}"

    def __len__(self):
        return len(self.design_seqs)

    def __getitem__(self, idx):
        return self.design_seqs[idx]

    def get_wt(self, seq:ProteinSequence, genename:str, read_from_file=False, filename=None):
        if read_from_file:
            seq = read_dna(os.path.join(self.dir, filename))
        self.wt = seq
        self.genename = genename

    def align(self, idx):
        design_seq = self.design_seqs[idx]
        alignment = align_optimal(design_seq, self.wildtype_seq)
        return alignment

    def get_design_sequence(self, seq, seq_id, score):
        self.design_seqs = seq
        self.design_seqs_id = seq_id
        self.design_seqs_score = score


    def get_msa(self, fasta_file):
        if not os.path.exists(os.path.join(self.dir, f"{self.genename}.msa")):
            subprocess.run(['jackhmmer', f"-A{os.path.join(self.dir, self.genename)}.msa",
                            fasta_file, UNIREF90])
        self.msa_file = f"{self.genename}.msa"

    def get_hmm(self):
        if not os.path.exists(os.path.join(self.dir, f"{self.genename}.hmm")):
            subprocess.run(['hmmbuild', "--hand", f"{os.path.join(self.dir, self.genename)}.hmm",
                            f"{os.path.join(self.dir,self.genename)}.msa"])
        self.hmm = HMM(f"{os.path.join(self.dir,self.genename)}.hmm")

    def get_pdb_sequence(self, pdbid):
        from biotite.database.rcsb import fetch
        self.pdbid = pdbid
        fetch(pdbid, 'mmcif', self.dir)
        fetch(pdbid, 'pdb', self.dir)
        self.pdb_file = os.path.join(self.dir, f"{pdbid}.pdb")
        self.mmcif_file = os.path.join(self.dir, f"{pdbid}.mmcif")


    def download_from_url(self, url, filename):
        if url not in self.url:
            self.url.append(url)
            self.url_filename.append(filename)
            response = requests.get(url, verify=False)
            # Check if the request was successful (status code 200)
            if response.status_code == 200:
                # Open a file in write-binary mode
                filename = f"{self.dir}/{filename}"
                with open(filename, 'wb') as file:
                    # Write the content of the response to the file
                    file.write(response.content)
                print(f"File downloaded successfully: {filename}")
            else:
                print(f"Failed to download file. Status code: {response.status_code}")

    def check_local_files(self, filename):
        if not os.path.exists(f"{self.dir}/{filename}"):
            return False
        return True