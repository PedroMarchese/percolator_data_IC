import os
import time
import threading as t
import numpy as np
import pandas as pd
from Bio import SeqIO

n_threads = 50
base_path = os.getcwd() or '/home/farminfo/PMarchese/percolator_data_IC'
checked = { 'found_dt': '', 'f_counter': 0, 'nf_counter': 0}
output_text = ''

class Protein:
    def __init__(self, identifier: str, sequence: str):
        self.identifier = identifier.replace('\n', '')
        self.sequence = sequence

    def get_string(self):
        return f'>{self.identifier}\n{self.sequence}\n'

def percentage_finder(annotated_p: pd.Series, start_sample: pd.Series, thread_id: int):
    a_counter = 0
    na_counter = 0
    chunk_output = ''
    nf_chunk_output = ''

    # ITERA O SUB ARRAY DE PROTEÍNAS ANOTADAS
    for protein in annotated_p:        
        # PROCURA NA SERIE O IDENTIFICADOR DA PROTEÍNA
        if start_sample.str.contains(protein.identifier).any():
            
            chunk_output += protein.get_string()
            a_counter += 1

            print(f'[T{thread_id}] | [{a_counter}/{annotated_p.size}] | {protein.identifier} encontrada!')
        else:
            nf_chunk_output += protein.get_string()
            na_counter += 1

            print(f'[T{thread_id}] | [{a_counter}/{annotated_p.size}] | {protein.identifier} não encontrada!')
            
    checked['f_counter'] += a_counter 
    checked['nf_counter'] += na_counter
    checked['found_dt'] += chunk_output
    checked['nfound_dt'] += nf_chunk_output
    
    print('='*15)
    print(f'{na_counter} NÃO ANOTADAS\n' + '-'*5)
    print(f'{a_counter / annotated_p.size} anotadas!')

def main():
    start_t = time.time()
    
    # LOAD FASTA ANNOTATED PROTEIN FILE    
    records = SeqIO.parse(f'{base_path}/data/mtb_proteome_cat.fasta', 'fasta')
    annotated_p = []
    # [annotated_p.append(Protein(rec.description, rec.seq)) for rec in records]
    for rec in records:
        if rec.description and rec.seq:
            annotated_p.append(Protein(rec.description, rec.seq))

    # LOAD ALL PROTEIN SAMPLE DATA (8 * 10^5)    
    with open(f'{base_path}/data/entire_output.txt', 'r') as f:
        lines = f.readlines()

    chunks = np.array_split(np.array(annotated_p), n_threads)    
    df = pd.DataFrame(chunks)
    #.replace('\n', '', regex=True)
    all_sample = pd.Series(lines)
    
    threads = []
    for i, df_chunk in df.iterrows():
        """
            CHUNKS = annotated protein list to splited array
            DF = all protein sample
        """
        threads.append(
            t.Thread(target=percentage_finder, args=(df_chunk, all_sample, i))
        )

    for thread in threads:
        thread.start()
    
    for thread in threads:
        thread.join()
        
    end_t = time.time()
    total_time = end_t - start_t
        
    # ESCRITA DOS OUTPUTS DE CADA THREAD
    with open(f'{base_path}/output/annotated.fasta', 'w', encoding='utf-8') as f:
        f.write(checked['found_dt'])

    with open(f'{base_path}/output/nannotated.fasta', 'w', encoding='utf-8') as f:
        f.write(checked['nfound_dt'])
        
    with open(f'{base_path}/output/report.txt', 'w', encoding='utf-8') as f:
        fc = checked['f_counter']
        nfc = checked['nf_counter']
        a_percentage = round(fc / len(annotated_p), 2)
        effectiveness = round(fc / all_sample.size, 2)
        
        report = f'[ENCONTRADAS] {checked}\n[NÃO ENCONTRADAS]{nfc}\n'
        report += f'[ANOTADAS_ENCONTRADAS] {a_percentage}%'
        report += f'[EFETIVIDADE_AMOSTRA] {effectiveness}%'
        report += f'[TEMPO_TOTAL] {total_time}\n[TEMPO/THREAD] {total_time/df.size}'
        f.write(report)
                
    print('Finished...')


main()