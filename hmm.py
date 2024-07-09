import biotite.database.entrez as entrez
import numpy as np
import subprocess
import pandas as pd
import os

from matplotlib import pyplot as plt
class HMM:
    def __init__(self, filename, description=None):
        self.filename = filename
        self.description = description
        self.trans_prob, self.match_prob, self.insert_prob = self.read_hmm_file()
        self.add_stop_codon()

    def subset(self, n):
        self.trans_prob = self.trans_prob.loc[:n, :]
        self.match_prob = self.match_prob.loc[:n, :]
        self.insert_prob = self.insert_prob.loc[:n, :]

    def add_stop_codon(self):
        self.match_prob.loc[self.match_prob.index[-1] + 1, :] = np.inf
        self.insert_prob.loc[self.insert_prob.index[-1] + 1, :] = self.insert_prob.loc[self.insert_prob.index[-1], :]
        self.trans_prob.loc[self.trans_prob.index[-1] + 1, :] = np.inf
        self.match_prob['*'] = np.inf
        self.match_prob.loc[self.match_prob.index[-1], '*'] = 0
        self.trans_prob.loc[self.trans_prob.index[-1], ['m->m', 'i->m', 'd->m', 'm->i', 'i->i']] = 0

    def read_hmm_file(self):
        start = False
        data = []
        with open(self.filename, 'r') as file:
            for line in file:
                if line.startswith('HMMER'):
                    continue
                elif line.startswith('HMM'):
                    start = True
                    data.append(line.strip().split())
                elif line.startswith('//') and start:
                    break
                elif start:
                    data.append(line.strip().split())

            trans_prob = []
            match_prob = []
            match_index = []
            ins_prob = []
            for i in range(len(data) // 3):
                trans_prob.append(data[3 * i + 1])
                match_index.append(data[3 * i + 2][0])
                match_prob.append(data[3 * i + 2][1:21])
                ins_prob.append(data[3 * i + 3])

            match_prob = pd.DataFrame(match_prob, index=match_index, columns=data[0][1:])
            trans_prob = pd.DataFrame(trans_prob[1:], columns=trans_prob[0],
                                     index=[x + 1 for x in range(len(trans_prob) - 1)])
            ins_prob = pd.DataFrame(ins_prob, index=match_index, columns=data[0][1:])
            trans_prob = trans_prob.apply(pd.to_numeric, args=('coerce',)).fillna(np.inf)
            match_prob = match_prob.apply(pd.to_numeric, args=('coerce',)).fillna(np.inf)
            ins_prob = ins_prob.apply(pd.to_numeric, args=('coerce',)).fillna(np.inf)
            ins_prob.reset_index(inplace=True)
            match_prob.reset_index(inplace=True)
            ins_prob.drop(columns='index', inplace=True)
            # row normalize the insertion probabilities to have std of 1 and mean of 0
            ins_mean = ins_prob.mean(axis=1)
            # match ins_mean dimensions to ins_prob
            ins_prob = ins_prob.apply(lambda x: x-x.mean(), axis=1)
            match_prob.drop(columns='index', inplace=True)
            match_prob.loc[0] = np.inf
            match_prob.loc[0,'M'] = 0
            return (trans_prob, match_prob, ins_prob)



