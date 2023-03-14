"""script to get the percentage of anotated smorfs in the percolator files"""

import os
import pandas as pd
import numpy as np
import threading
import time as t
from time import time
from random import random

data_filename = 'all_results_norepeats.txt'
n_threads = 30

def protein_gather(array: np.array, thread_id, base_path):
    start = time()
    output_ids = []
    
    counter = 0
    for proteins in array:
        if counter % 5000 == 0:
            print(f'[THREAD_ID] {thread_id} | [COUNTER] {counter}/{array.size} | {round(counter/array.size*100, 0)}\n')

        protein_array = proteins.split(',')
        output_ids = output_ids + protein_array

        counter += 1
        
    with open(f'{base_path}/proteinIds_outputs/{thread_id}.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(output_ids))
        
    with open(f'{base_path}/entire_output.txt', 'a', encoding='utf-8') as f:
        f.write('\n'.join(output_ids) + '\n')

    print(f'[END] Thread {thread_id} finished in {round(time()-start, 0)} seconds!')


def main():
    start = time()
    
    # base_path = '/home/farminfo/ACanedo/percolator/'
    base_path = ''
    if not base_path:
        base_path = os.getcwd()

    # Threads variables
        threads = []
            
    #Reading file
    df = pd.read_csv(f'{base_path}/{data_filename}', sep="\t")
    np_array = np.array(df['proteinIds'].values)
    chunks = np.array_split(np_array, n_threads)
        
    for i, chunk in enumerate(chunks):
        threads.append(threading.Thread(target=protein_gather, args=(chunk, i, base_path)))
        
    for t in threads:
        t.start()
        
    for t in threads:
        t.join()

    print(f'{time()-start} segundos de execução. Isso é tudo, pessoal!')
    
    
main()