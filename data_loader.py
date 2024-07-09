from biotite.sequence import ProteinSequence
from biotite.sequence.align import align_optimal
import subprocess
import requests
import os
from hmm import HMM
from constants import UNIREF50
# define dataset class
class Dataset:
    def __init__(self, name):
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

    def get_wt(self, seq:ProteinSequence, genename:str):
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
                            fasta_file, UNIREF50])
        self.msa_file = f"{self.genename}.msa"

    def get_hmm(self):
        if not os.path.exists(os.path.join(self.dir, f"{self.genename}.hmm")):
            subprocess.run(['hmmbuild', "--hand", f"{os.path.join(self.dir, self.genename)}.hmm",
                            f"{os.path.join(self.dir,self.genename)}.msa"])
        self.hmm = HMM(f"{os.path.join(self.dir,self.genename)}.hmm")


    def get_protein_sequence(self, idx):
        pass

    def get_wildtype_sequence(self):
        pass

    def get_pdb_sequence(self):
        pass

    def download_from_url(self, url, filename):
        if url not in self.url:
            self.url.append(url)
            self.url_filename.append(filename)
            response = requests.get(url)
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
