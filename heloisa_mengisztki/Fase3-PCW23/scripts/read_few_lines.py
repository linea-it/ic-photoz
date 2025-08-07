import argparse
import pyarrow.parquet as pq
import os

BUFFER_SIZE = 1<<20

def read_write_parquet(input_path, output_path, rows):
    batch_size = 10  

    with pq.ParquetFile(input_path, buffer_size=BUFFER_SIZE) as file:
        batches = file.iter_batches(batch_size=rows, use_pandas_metadata=True)

        block = next(iter(batches))
        block.to_pandas().to_parquet(output_path)
  
    print("===> arquivo salvo: ", output_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Ler e gravar arquivo Parquet')
    parser.add_argument('input_path', type=str, help='Caminho para o arquivo Parquet de entrada')
    parser.add_argument('output_path', type=str, help='Caminho para o arquivo Parquet de saída')
    parser.add_argument('rows', type=int, help='Número de linhas a serem lidas')
    
    args = parser.parse_args()
    
    read_write_parquet(args.input_path, args.output_path, args.rows)
