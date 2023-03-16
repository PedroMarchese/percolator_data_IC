import os
import pandas as pd
import numpy as np
import threading as t
from Bio import SeqIO

n_threads = 20
base_path = os.getcwd() or '/home/farminfo/PMarchese/percolator_data_IC'
checked = { 'data': '', 'counter': 0}

class Protein:
    def __init__(self, identifier: str, sequence: str):
        self.identifier = identifier.replace('\n', '')
        self.sequence = sequence

    def get_string(self):
        return f'>{self.identifier}\n{self.sequence}\n'

def percentage_finder(all_proteins: pd.Series, thread_id: int):
    
    records = SeqIO.parse(f'{base_path}/data/mtb_proteome_cat.fasta', 'fasta')
    found_proteins = []

    annotated_counter = 0
    for i, record in enumerate(records):
        identifiers = str(record.description)
        sequences = str(record.seq)

        output_string = ''
        for protein in found_proteins:
            output_string += protein.get_string()
            print(output_string)
        
        # TRATA ESSES DADOS
        if all_proteins.str.contains(identifiers).any():
            print(f'[{annotated_counter} / {all_proteins.size}] | [T{thread_id}] | {identifiers} encontrado!')
            
            found_proteins.append(Protein(identifiers, sequences))
            annotated_counter += 1
            
    checked['counter'] += annotated_counter 
    checked['data'] += output_string
    
    print(f'{round(annotated_counter / all_proteins.size * 100, 3)}% anotadas!')

def main():
    # LOAD FASTA ANNOTATED PROTEIN FILE    
    records = SeqIO.parse(f'{base_path}/data/mtb_proteome_cat.fasta', 'fasta')
    
    with open(f'{base_path}/data/entire_output.txt', 'r') as f:
        lines = f.readlines()

    annotated = pd.Series(lines).replace('\n', '', regex=True)
    for rec in records:
        print(rec)
    # chunks = np.array_split(np.array(list(records)), n_threads)
    # print(chunks)
    # df = pd.DataFrame(chunks).replace('\n', '', regex=True)

    # print(df)

    # threads = []
    # for i, row in df.iterrows():
    #     threads.append(
    #         t.Thread(target=percentage_finder, args=(row, i))
    #     )

    # for thread in threads:
    #     thread.start()
    
    # for thread in threads:
    #     thread.join()
        
    # print('Finished...')


main()