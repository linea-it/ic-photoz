#!/usr/bin/env python3

import os
import glob
from sys import argv
from pyarrow.parquet import ParquetFile
import numpy as np

dp0_dir = '/lustre/t1/cl/lsst/dp0.2/'
helo = '/lustre/t0/scratch/users/heloisa.mengisztki/'

dir50  = helo+'dp0.2_pre_processed_50_original_files'
dir100 = helo+'dp0.2_pre_processed_100_original_files'
dir150 = helo+'dp0.2_pre_processed_150_original_files'

print(50) 
os.system(f'du -hs {dir50}')
print(100) 
os.system(f'du -hs {dir100}')
print(150) 
os.system(f'du -hs {dir150}')

