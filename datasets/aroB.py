import pandas as pd
from data_loader import Dataset
from biotite.sequence import ProteinSequence, NucleotideSequence
import numpy as np
import os

aroB = Dataset('aroB')
for f in ['aroB.fa', 'infA_aroB_20240411.csv']:
    aroB.check_local_files(f)

aroB.get_wt('', 'aroB', read_from_file=True, filename='aroB.fa')
aroB.get_msa(os.path.join(aroB.dir, 'aroB.fa'))
aroB.get_hmm()

aroB.get_pdb_sequence('5eks')

def aroB_enrichment(design_seqs, subset):
    if subset in ['ERP017','ERP766']:
        design_seqs = design_seqs[design_seqs['aroB_location']==subset]
    elif subset == 'all':
        pass
    else:
        raise ValueError('subset must be ERP017, ERP766 or all')
    full_seq = design_seqs['full_seq'].values
    prot_seq = [NucleotideSequence(seq).translate(complete=True) for seq in full_seq]
    variant_id = design_seqs['unique_seq_id'].values
    # check that the log enrichement is the log of the enrichment
    assert (np.log(design_seqs['enrichment_aroB_1'][:50])-design_seqs['log_enrichment_aroB_1'][:50].astype(float) < 1e-6).sum() == 50
    assert (np.log(design_seqs['enrichment_aroB_2'][:50])-design_seqs['log_enrichment_aroB_2'][:50].astype(float) < 1e-6).sum() == 50
    assert (np.log(design_seqs['enrichment_aroB_3'][:50])-design_seqs['log_enrichment_aroB_3'][:50].astype(float) < 1e-6).sum() == 50
    for i in [1,2,3]:
        design_seqs.loc[:,f'log_enrichment_aroB_{i}'] = np.log(design_seqs[f'enrichment_aroB_{i}']+1e-20)

    # not going to care about the wt enrichment for now
    # plot the correlation between the log enrichement for the 3 assays
    """
    The following section checks the correlation between the log enrichments for the 3 assays
    The correlation is high between the assays and we use the average of the 3 assays as the target
    """
    # import seaborn as sns
    # import matplotlib.pyplot as plt
    # # make a figure with 3 subplots
    # fig, ax = plt.subplots(1,3, figsize=(20, 8))
    # # plot the scatter plots
    # sns.scatterplot(x='log_enrichment_aroB_1', y='log_enrichment_aroB_2', data=design_seqs, ax=ax[0])
    # sns.scatterplot(x='log_enrichment_aroB_1', y='log_enrichment_aroB_3', data=design_seqs, ax=ax[1])
    # sns.scatterplot(x='log_enrichment_aroB_2', y='log_enrichment_aroB_3', data=design_seqs, ax=ax[2])
    # plt.show()
    # plt.close()
    scores = design_seqs[[f'log_enrichment_aroB_{i}' for i in [1,2,3]]].mean(axis=1)
    return prot_seq, variant_id, scores

design_seqs = pd.read_csv(os.path.join(aroB.dir, 'infA_aroB_20240411.csv'))
# wt_info = design_seqs[design_seqs['unique_seq_id']=='seq_0001']
# wt_info['enrichment_aroB_3']
# design_seqs.columns
"""
drop sequences that do not contain entangled sequences 
"""
design_seqs = design_seqs.dropna(subset=['full_seq'])
"""
add design sequence to object 
"""
aroB.get_design_sequence(*aroB_enrichment(design_seqs, 'all'))
# there are two subsets of sequences in the csv file based on aroB_location
aroB.design_seqs_score