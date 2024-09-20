from data_loader import Dataset
import os
import pandas as pd
import numpy as np
from biotite.sequence import ProteinSequence, NucleotideSequence


infA = Dataset('infA')
for f in ['infA_aroB_20240411.csv','infA.fa']:
    infA.check_local_files(f)

infA.get_wt('', 'infA', read_from_file=True, filename='infA.fa')
infA.get_msa(os.path.join(infA.dir, 'infA.fa'))
infA.get_hmm()

infA.get_pdb_sequence('2n78')

def infA_enrichment(design_seqs, subset):
    if subset in ['ERP017','ERP766']:
        design_seqs = design_seqs[design_seqs['aroB_location']==subset]
    elif subset == 'all':
        pass
    else:
        raise ValueError('subset must be ERP017, ERP766 or all')
    full_seq = design_seqs['full_seq'].values
    # extract the capitalized letters from the full sequence
    entangled_seq = ["".join([x for x in seq if x.isupper()]) for seq in full_seq]
    entangled_seq = [x[1:] for x in entangled_seq]
    prot_seq = [NucleotideSequence(seq).translate(complete=True)[:-1] for seq in entangled_seq]
    variant_id = design_seqs['unique_seq_id'].values
    # check that the log enrichement is the log of the enrichment
    assert (np.log(design_seqs['enrichment_infA_1'][:20])-design_seqs['log_enrichment_infA_1'][:20].astype(float) < 1e-6).sum() == 20
    assert (np.log(design_seqs['enrichment_infA_2'][:20])-design_seqs['log_enrichment_infA_2'][:20].astype(float) < 1e-6).sum() == 20
    assert (np.log(design_seqs['enrichment_infA_3'][:20])-design_seqs['log_enrichment_infA_3'][:20].astype(float) < 1e-6).sum() == 20
    for i in [1,2,3]:
        design_seqs.loc[:,f'log_enrichment_infA_{i}'] = np.log(design_seqs[f'enrichment_infA_{i}']+1e-20)

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
    # sns.scatterplot(x='log_enrichment_infA_1', y='log_enrichment_infA_2', data=design_seqs, ax=ax[0])
    # sns.scatterplot(x='log_enrichment_infA_1', y='log_enrichment_infA_3', data=design_seqs, ax=ax[1])
    # sns.scatterplot(x='log_enrichment_infA_2', y='log_enrichment_infA_3', data=design_seqs, ax=ax[2])
    # plt.show()
    # plt.close()
    scores = design_seqs[[f'log_enrichment_infA_{i}' for i in [1,2,3]]].mean(axis=1)
    return prot_seq, variant_id, scores

design_seqs = pd.read_csv(os.path.join(infA.dir, 'infA_aroB_20240411.csv'))
wt_info = design_seqs[design_seqs['unique_seq_id']=='seq_0001']
wt_enrichment = wt_info[['enrichment_infA_1', 'enrichment_infA_2', 'enrichment_infA_3', 'log_enrichment_infA_1', 'log_enrichment_infA_2', 'log_enrichment_infA_3']]
design_seqs = design_seqs.dropna(subset=['full_seq'])
"""
drop sequences that do not contain entangled sequences 
"""
design_seqs = design_seqs.dropna(subset=['full_seq'])
"""
add design sequence to object 
"""
infA.get_design_sequence(*infA_enrichment(design_seqs, 'all'))
# there are two subsets of sequences in the csv file based on aroB_location
