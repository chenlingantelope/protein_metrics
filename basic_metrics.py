# metrics that are independent from external information
# sequence length, aa composition, charge, hydrophobicity, aromaticity, polarity
import pandas as pd


def sequence_length(seq):
    return len(seq)

def aa_composition(seq):
    aa = 'ACDEFGHIKLMNPQRSTVWY-'
    aa_count = {i: 0 for i in aa}
    for x in seq:
        aa_count[x.upper()] += 1
    return pd.Series(aa_count)

def hamming_distance(seq, ref_seq):
    seq = ''.join([x for x in seq if x.isupper() or x == '-'])
    ref_seq = ''.join([x for x in ref_seq if x.isupper() or x == '-'])
    return sum([i != j for i, j in zip(seq, ref_seq)])


def sequence_charge(seq, pH):
    pka = {
        'AA': ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
        'N_pKa': [9.69, 9.04, 8.80, 9.60, 10.28, 9.67, 9.13, 9.60, 9.17, 9.68,
                  9.60, 8.95, 9.21, 9.13, 10.96, 9.15, 9.62, 9.39, 9.11, 9.62],
        'C_pKa': [2.34, 2.17, 2.02, 1.88, 1.96, 2.19, 2.17, 2.34, 1.82, 2.36,
                  2.36, 2.18, 2.28, 2.20, 1.99, 2.21, 2.11, 2.38, 2.20, 2.32],
        'R_pKa': [None, 12.48, None, 3.65, 8.18, 4.25, None, None, 6.00, None,
                  None, 10.53, None, None, None, None, None, None, 10.07, None]
    }
    pka = pd.DataFrame(pka)
    pka.set_index('AA', inplace=True)
    # convert pka to a dictionary using AA as key and N_pKa, C_pKa, R_pKa as values
    pka_dict = pka.to_dict(orient='index')
    total_charge = 0
    h = 10**(-pH)
    # remove gaps and convert to uppercase
    seq = ''.join([x.upper() for x in seq if x!='-'])
    for i, x in enumerate(seq):
        charge = 0
        if not pd.isna(pka_dict[x]['R_pKa']):
            if x in 'HKR':
                charge += h / (h + 10**(-pka_dict[x]['R_pKa']))
            elif x in 'DECY':
                charge += h / (h + 10**(-pka_dict[x]['R_pKa'])) - 1
        if i == 0:
            charge += h / (h + 10**(-pka_dict[x]['N_pKa']))
        elif i == len(seq) - 1:
            charge += h / (h + 10**(-pka_dict[x]['C_pKa'])) - 1
        total_charge += charge
    return total_charge


def aa_properties(seq):
    seq = ''.join([x.upper() for x in seq if x!='-'])
    n_positively_charged = sum([x in 'HKR' for x in seq])
    n_negatively_charged = sum([x in 'DE' for x in seq])
    n_polar = sum([x in 'STNQCY' for x in seq])
    n_hydrophobic = sum([x in 'GAVLIPFMW' for x in seq])
    n_aromatic = sum([x in 'HFYW' for x in seq])
    return pd.Series([n_positively_charged, n_negatively_charged, n_polar, n_hydrophobic, n_aromatic],
                     index=['n_positively_charged', 'n_negatively_charged', 'n_polar', 'n_hydrophobic', 'n_aromatic'])


def all_basic_metrics(design_seqs, wildtype_seq, pH = 7.4):
    # define pandas dataframe for basic metrics, indexed by sequence id in design_seqs
    basic_metrics = pd.DataFrame(index=[seq.id for seq in design_seqs])
    # sequence length
    basic_metrics['len'] = [sequence_length(seq.seq) for seq in design_seqs]
    # amino acid composition
    aa_counts = pd.concat([aa_composition(seq.seq) for seq in design_seqs], axis=1).T
    basic_metrics[aa_counts.columns] = aa_counts.values
    # compute hamming distance to wildtype (within the position of the gene profile)
    basic_metrics['hamming_distance'] = [hamming_distance(seq, wildtype_seq) for seq in design_seqs]
    # sequence charge
    # https://www.bu.edu/aldolase/biochemistry/html_docs/Calculating_pI.ppt#:~:text=Multiply%20the%20proportion%20charged%20by,to%20get%20the%20net%20charge.
    basic_metrics['charge'] = [sequence_charge(seq, pH) for seq in design_seqs]
    # amino acid properties
    aa_prop = pd.concat([aa_properties(seq.seq) for seq in design_seqs], axis=1).T
    basic_metrics[aa_prop.columns] = aa_prop.values
    return basic_metrics