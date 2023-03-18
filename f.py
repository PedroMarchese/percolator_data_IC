import os
import pandas as pd
import numpy as np
import threading as t
from Bio import SeqIO

n_threads = 10
base_path = os.getcwd() or '/home/farminfo/PMarchese/percolator_data_IC'

class Protein:
    def __init__(self, identifier: str, sequence: str):
        self.identifier = identifier.replace('\n', '')
        self.sequence = sequence

    def get_string(self):
        return f'>{self.identifier}\n{self.sequence}\n'

def percentage_finder(all_proteins: pd.Series):
    records = SeqIO.parse(f'{base_path}/data/mtb_proteome_cat.fasta', 'fasta')
    found_proteins = []


def main():
        
    df = pd.read_csv(f'{base_path}/data/entire_output.txt', header=0)
    anno = df[df['proteinId'].str.contains('ANNO')]
    mask = anno.duplicated(keep=False)
    unique_anno = anno[~mask]

    unique_anno.to_csv(f'{base_path}/output/filtered_anno.csv', index=False)

main()