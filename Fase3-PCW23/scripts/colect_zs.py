import argparse
import pyarrow.parquet as pq
import os
import pandas as pd
import h5py
import sys

def read_hdf5_file(file_path, zmodes):       
    with h5py.File(file_path, 'r') as file:
        sys.stdout.write(".")
        sys.stdout.flush()
        
        zmode_mean = file['ancil']['zmode'][:].mean()
        #zmode_mean = file['ancil']['zmean'][:]
        zmodes.append(zmode_mean)

def main(output_path, base_path, test_id):
    print(f"processing - {test_id}")
        
    folder_path = f"{base_path}/{test_id}/output/"
    for root, _, files in os.walk(folder_path):
        zmodes = []
        files_names = []
        
        for file_name in files:
            if file_name.endswith('.pq'):
                file_path = os.path.join(root, file_name)
                files_names.append(file_name)
                read_hdf5_file(file_path,zmodes)
                
        data = {
            'zmode_mean': zmodes,
            'filepath': files_names
        }
    
        tasks_df = pd.DataFrame(data=data)
        output_file = f"{output_path}/{test_id}_output.csv"
        tasks_df.to_csv(output_file, index=False, header=['zmode_mean', 'filepath'], mode='w')
        
        print("Qtd arquivos - ", len(zmodes))
                

if __name__ == "__main__":
    output_path = "~/ic-photoz/Fase3-PCW23/outputs"
    base_path = "/lustre/t1/cl/lsst/pz_project/pcw_2023" #/lustre/t0/scratch/users/julia/pcw_2023
    tests = [
        #t1
        "fzboost_trunc4_chunk_150k",
        "fzboost_all_dec_cases_chunk_150k"
    ]
    
    for test_id in tests:
        main(output_path, base_path, test_id)