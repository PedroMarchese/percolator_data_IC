import os
import pandas as pd
import numpy as np
import threading as t
from Bio import SeqIO

n_threads = 30
base_path = os.getcwd() or '/home/farminfo/PMarchese/percolator_data_IC'
checked = { 'data': '', 'counter': 0}
output_text = ''

class Protein:
    def __init__(self, identifier: str, sequence: str):
        self.identifier = identifier.replace('\n', '')
        self.sequence = sequence

    def get_string(self):
        return f'>{self.identifier}\n{self.sequence}\n'

def percentage_finder(proteins: pd.Series, all_proteins: pd.Series, thread_id: int):
    found_proteins = []

    annotated_counter = 0
    for i, protein in enumerate(proteins):
        output_string = ''
        for protein in found_proteins:
            output_string += protein.get_string()
            print(output_string)
        
        # TRATA ESSES DADOS
        # if all_proteins.str.contains(protein.description).any():
        return print(all_proteins)
        if all_proteins.contains(protein.description):
            print(f'[{annotated_counter} / {all_proteins.size}] | [T{thread_id}] | {protein.identifier} encontrado!')
            
            output_text += protein.get_string()
            annotated_counter += 1
            
    checked['counter'] += annotated_counter 
    checked['data'] += output_string
    
    print(f'{round(annotated_counter / all_proteins.size * 100, 3)}% anotadas!')

def main():
    # LOAD FASTA ANNOTATED PROTEIN FILE    
    records = SeqIO.parse(f'{base_path}/data/mtb_proteome_cat.fasta', 'fasta')
    proteins = []
    [proteins.append(Protein(rec.description, rec.seq)) for rec in records]
    
    with open(f'{base_path}/data/entire_output.txt', 'r') as f:
        lines = f.readlines()

    chunks = np.array_split(np.array(proteins), n_threads)
    df = pd.DataFrame(chunks).replace('\n', '', regex=True)
    annotated = pd.Series(lines).replace('\n', '', regex=True)

    

    threads = []
    for i, chunk in df.iterrows():
        threads.append(
            t.Thread(target=percentage_finder, args=(, df, i))
        )
        break

    for thread in threads:
        thread.start()
    
    for thread in threads:
        thread.join()
        
    print('Finished...')


main()