import argparse
import os
import pandas as pd
import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt

column_names = [f"bin_{round(value, 2)}" for value in np.linspace(0,3,300)]

def get_host(file_name, test_id):
    
    file_path = f"~/ic-photoz/Fase3-PCW23/results/tests/{test_id}.csv"
    
    df = pd.read_csv(file_path)
    condition = df['file_name'].str.contains(file_name.replace(".pq", "")) 
    filtered_df = df[condition]
    host = filtered_df['host'].values
    return host[0]

def read_hdf5_file(output_file, file_path, file_name, test_id):       
    with h5py.File(file_path, 'r') as file:
        sys.stdout.write(".")
        sys.stdout.flush()
        
        zmode = file['ancil']['zmode'][:]
        bins = file['meta']['xvals'][:][0]
        
        n, bins, patches = plt.hist(zmode, bins=bins)
        
        #print(len(n), len(bins))
        
        df_row = pd.DataFrame([n])
        df_row.loc[0, 'host']=get_host(file_name, test_id)
        df_row.loc[0, 'file_name']=file_name
        df_row.to_csv(output_file, index=False, header=False, mode='a')

def main(output_path, base_path, test_id):
    print(f"processing - {test_id}")
        
    folder_path = f"{base_path}/{test_id}/output/"
    for root, _, files in os.walk(folder_path):
        
        column_names.append('host')
        column_names.append('file_name')        
        
        output_file = f"{output_path}/nz_zmodes_hist_{test_id}.csv"
        df = pd.DataFrame(columns=column_names)
        df.to_csv(output_file, index=False, header=True, mode='w')
        
        for file_name in files:
            if file_name.endswith('.pq'):
                file_path = os.path.join(root, file_name)
                read_hdf5_file(output_file, file_path, file_name, test_id)
                  

if __name__ == "__main__":
    output_path = "~/ic-photoz/Fase3-PCW23/outputs"
    base_path = "/lustre/t1/cl/lsst/pz_project/pcw_2023" #/lustre/t0/scratch/users/julia/pcw_2023
    tests = [
        #t1
        #"fzboost_trunc4_chunk_150k",
        "fzboost_all_dec_cases_chunk_150k_100x",#nop
        "fzboost_trunc4_chunk_150k_100x",#nop
        "bpz_trunc4_chunk_150k",
        "bpz_all_dec_cases_chunk_150k",
        "bpz_seds_default",
        "bpz_seds_2x_default",
        "test_hardware_t1"
    ]
    
    for test_id in tests:
        main(output_path, base_path, test_id)