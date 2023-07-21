import os
import yaml
import csv
import argparse

from pyarrow.parquet import ParquetFile

def process_files(folder_path, output_path):       
    for root, directories, files in os.walk(folder_path):
        preprocess_dirs = [dir_name for dir_name in directories if "pre_processed" in dir_name.lower()]

        for preprocess_dir in preprocess_dirs:
            dir_path = os.path.join(root, preprocess_dir)  
            output_file = os.path.join(output_path, f"{preprocess_dir}_data.csv") 
            
            print("Processing -", dir_path)
            
            with open(output_file, 'a', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)

                if os.path.getsize(output_file) == 0:
                    csv_writer.writerow(['dir', 'file_name', 'row_counts', 'file_size'])
                
                i = 0
                
                for file_name in os.listdir(dir_path):
                    print(i + 1)
                    if file_name.endswith('.pq'):
                        file_path = os.path.join(dir_path, file_name)
                        size = os.path.getsize(file_path)
                        rows = 0
                        with ParquetFile(file_path) as f:
                            rows = f.metadata.num_rows
                        csv_writer.writerow([file_path, file_name, rows, size])
                    i += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_path", help="Path to the folder to scan for YAML files.")
    parser.add_argument("output_path", help="Path to the output CSV file.")
    args = parser.parse_args()
    
    process_files(args.folder_path, args.output_path)