import os
import pandas as pd
from data_loader import Dataset
from biotite.sequence import ProteinSequence, NucleotideSequence
from biotite.sequence.io import fasta
from biotite.application.blast import BlastWebApp
from biotite.database import entrez
from tempfile import gettempdir


def mutate_wt(wt, x):
    if x == '_wt':
        return wt
    else:
        x = x.split('(')[1].split(')')[0]
    if 'ins' in x:
        pos = int(x.split('_')[0][1:])
        ins = x.split('ins')[1]
        assert wt[pos-1] == x[0]
        return wt[:pos-1] + ProteinSequence(ins) + wt[pos-1:]
    elif 'del' in x:
        x = x.replace('del', '')
        if '_' in x:
            x1 = int(x.split('_')[0][1:])
            x2 = int(x.split('_')[1][1:])
            assert wt[x1-1] == x.split('_')[0][0]
            assert wt[x2-1] == x.split('_')[1][0]
            return wt[:x1-1] + wt[x2-1:]
        else:
            x1 = int(x[1:])
            assert wt[x1-1] == x[0]
            return wt[:x1-1] + wt[x1:]
    else:
        pos = int(x[1:-1])
        assert wt[pos-1] == x[0]
        return wt[:pos-1] + ProteinSequence(x[-1]) + wt[pos:]


dimple = Dataset('dimple')
urls = ['https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-023-02880-6/MediaObjects/13059_2023_2880_MOESM1_ESM.csv',
        'https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-023-02880-6/MediaObjects/13059_2023_2880_MOESM2_ESM.txt',
        'https://raw.githubusercontent.com/odcambc/DIMPLE/master/tests/Kir.fa']

url_files = ['seq.csv', 'score.txt', 'wt.fa']
for i,x in enumerate(urls):
    dimple.download_from_url(x, url_files[i])

# read the wt.fa file
seq = fasta.FastaFile.read(os.path.join(dimple.dir, url_files[2]))
seq = fasta.get_sequence(seq)
wt = seq[582:1893].translate(complete=True)
dimple.get_wt(wt, 'Kir21')

# mutate the wildtype sequence
data = pd.read_csv(os.path.join(dimple.dir, url_files[1]), sep='\t')
variants = [mutate_wt(wt, x) for x in data['hgvs']]
variantid = list(data['hgvs'])
dimple.get_design_sequence(variants, variantid, data['score'])

# get msa alignment
dimple.get_msa(os.path.join(dimple.dir, 'wt.fa'))
dimple.get_hmm()
