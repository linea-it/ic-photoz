from collections import namedtuple
from datetime import datetime
from re import search
from sys import argv, stdin
import yaml
import pandas as pd
import os

import sys
from itertools import islice

JOB_EXECUTING = \
        r'[.](\d+)[.].*(\d{2}/\d{2} \d{2}:\d{2}:\d{2}) Job executing.*<(.*):'

JOB_TERMINATED = \
    r'[.](\d+)[.].*(\d{2}/\d{2} \d{2}:\d{2}:\d{2}) Job terminated'

CHUNKS = \
    r"on chunk 0 - (\d+)"

OUTPUT_CHUNK_LOG = \
    r"output_estimate: output/(.*?)\.pq"

TaskData = namedtuple('TaskData', ('task_id', 'time_begin', 'time_end', 'time_diff', 'file_name', 'chunks', 'size'))

def process_output_file(task_id, path_logs_files):
    for root, _, files in os.walk(path_logs_files):
        for file_name in files:
            if file_name.endswith(f'-{task_id}.out'):
                sys.stdout.write(".")
                sys.stdout.flush()
                
                file_path = os.path.join(path_logs_files, file_name)
                
                file_name = ""
                chunks = 0
                file_size = 0
                
                with open(file_path) as f:
                    for line in f:
                        m = search(CHUNKS, line)
                        if m:
                            chunks =int(m.group(1))
                        
                        m = search(OUTPUT_CHUNK_LOG, line)
                        if m:
                            file_name = m.group(1)[len("inprogress_"):]
                            print(file_path)
                            file_size = os.path.getsize(file_path)
                          
                return file_name, chunks, file_size



def process_file(file_path, root, process_id):
    print(f"Processando o arquivo: {file_path}")
    
    path_logs_files = os.path.join(root, "log/")
    
    task_map = {}
    host_map = {}
    
    f = open(file_path)
    
    for line in f:
        m = search(JOB_EXECUTING, line)
        if m:
            task_id = int(m.group(1))
            time_begin = datetime.strptime(m.group(2), '%m/%d %X')
            host = m.group(3)
            task_map[task_id] = (time_begin, host)
            
        m = search(JOB_TERMINATED, line)
        if m:
            task_id = int(m.group(1))
            time_end = datetime.strptime(m.group(2), '%m/%d %X')
            time_begin = task_map[task_id][0]
            host = task_map[task_id][1]
            time_diff = (time_end - time_begin).total_seconds()
            host_map.setdefault(host, [])
            
            file_name, chunks, file_size = process_output_file(task_id, path_logs_files)
            
            t = TaskData(task_id, time_begin, time_end, time_diff, file_name, chunks, file_size)
            host_map[host].append(t)
            
    host_list = [(k, v) for (k, v) in host_map.items()]

    csv_file_path = f'~/{process_id}.csv'
        
    df = pd.DataFrame(columns=['task_id', 'time_begin', 'time_end', 'time_diff', 'file_name', 'chunks', 'size', 'host'])
    df.to_csv(csv_file_path, index=False, mode='a')
    
    for host, tasks in host_list:
        tasks_df = pd.DataFrame(tasks)
        tasks_df['host'] = host
        tasks_df.to_csv(csv_file_path, index=False, header=False, mode='a')

def process_files_in_folder(folder_path, process_id):
    for root, _, files in os.walk(folder_path):
        for file_name in files:
            if file_name.endswith('rail-condor.log'):
                file_path = os.path.join(root, file_name)
                process_file(file_path, root, process_id)
                

def main(csv_file):
    
    df = pd.read_csv(csv_file)
    folders = df['file_path'].tolist()

    for row in df[['file_path', 'process_id']].itertuples(index=False):
        folder_path, process_id = row
        if os.path.exists(folder_path):
            process_files_in_folder(folder_path, process_id)
        else:
            print(f"A pasta '{folder_path}' não existe.")

def process_for_folder():
    folder_path = "/home/henrique.almeida/results/"
    process_id = "log-10B-flexzboost-2" #log-10B-flexzboost/   log-10B-flexzboost-2/
    if os.path.exists(folder_path):
        process_files_in_folder(folder_path, process_id)
    else:
        print(f"A pasta '{folder_path}' não existe.")
            
if __name__ == "__main__":
    csv_file = "~/ic-photoz/Fase3-PCW23/results/results_1.csv"
    main(csv_file)
    #process_for_folder()