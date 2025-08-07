#!/usr/bin/env python3

import os
import glob
from sys import argv
from pyarrow.parquet import ParquetFile
import numpy as np

dp0_dir = '/lustre/t1/cl/lsst/dp0.2/'
helo = '/lustre/t0/scratch/users/heloisa.mengisztki/'

dir50  = helo+'dp0.2_pre_processed_50_original_files/'
dir100 = helo+'dp0.2_pre_processed_100_original_files/'
dir150 = helo+'dp0.2_pre_processed_150_original_files/'

files50  = glob.glob(dir50+'/*')
files100 = glob.glob(dir100+'/*')
files150 = glob.glob(dir150+'/*')

all_files = [files50, files100, files150] 

'''
for files in all_files: 
    counts = []
    for i in files: 
        sum_rows = 0
        with ParquetFile(i) as f:
            #print(i, f.metadata.num_rows)
            sum_rows += f.metadata.num_rows
        #print(d, sum_rows)
        counts.append(sum_rows)
    total = sum(counts) 
    print(f'total: {total} rows')

'''
print()
dp0_files = glob.glob(dp0_dir+'*.parq') 
counts = []
for i in dp0_files: 
    sum_rows = 0
    with ParquetFile(i) as f:
        #print(i, f.metadata.num_rows)
        sum_rows += f.metadata.num_rows
    #print(d, sum_rows)
    counts.append(sum_rows)
total = sum(counts) 
print(f'total dp0 original: {total} rows')





