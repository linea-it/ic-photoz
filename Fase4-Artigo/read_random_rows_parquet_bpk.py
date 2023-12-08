import pandas as pd
import numpy as np

from os.path import basename, join, splitext
from pyarrow.parquet import ParquetFile
from pyarrow import ArrowInvalid
from sys import argv

BUFFER_SIZE = 1<<20
DEFAULT_EXT = '.csv'
TEMPLATE = "sample_{fname}{ext}"

def build_output_name(file_name, template, extension=DEFAULT_EXT):
    file_name = basename(file_name)
    split = splitext(file_name)
    ext = extension or split[1]

    return template.format(fname=split[0],  ext=ext)

def read_random_file(parquet_file, output_dir, num_rows):

   try:
        with ParquetFile(parquet_file, buffer_size=BUFFER_SIZE) as f:
            parquet_data = pd.read_parquet(parquet_file)
            sample_data = parquet_data.sample(n=num_rows)

            output_name = build_output_name(parquet_file, TEMPLATE) 
            output_path = join(output_dir, output_name)

            sample_data.to_csv(output_path, index=False)

   except ArrowInvalid as e:
        print('Error: parquet library returned error.', file=stderr)
        print(e, file=stderr)
        raise SystemExit(1)
   except OSError as e:
        print('Error: %s' % e, file=stderr)
        raise SystemExit(1)


def main():
    parquet_file = argv[1]
    output_dir = argv[2]
    num_rows = int(argv[3])

    read_random_file(parquet_file, output_dir, num_rows)

if __name__ == '__main__': main()