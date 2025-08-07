import os
import yaml
import csv
import argparse

## python script_1.py caminho_da_pasta ~/output.csv

def open_yaml_file(file_path, csv_writer, folder_path):
    try:
        with open(file_path, 'r') as file:
            yaml_data = yaml.safe_load(file)
            
            print(f"Reading file - {file_path}")
            
            comment = yaml_data.get('comment')
            process_id = yaml_data.get('process_id')
            time = yaml_data.get('time stats')
            start = time.get('Process started')
            end = time.get('Process finished')
            duration = time.get('Total duration')
            user = yaml_data.get('user')            
            
            csv_writer.writerow([comment, process_id, start, end, duration, user, folder_path])
    except Exception as e:
        print(f"Erro ao processar o arquivo {file_path}: {e}")

        
def process_yaml_files(folder_path, output_file):
        
    csv_exists = os.path.isfile(output_file)

    with open(output_file, 'a', newline='') as csv_file: #use a to append
        csv_writer = csv.writer(csv_file)

        if not csv_exists:
            csv_writer.writerow(['comment', 'process_id', 'start', 'end', 'duration', 'user', 'file_path'])

        for root, _, files in os.walk(folder_path):
            for file_name in files:
                if file_name.endswith('process_info.yaml'):
                    file_path = os.path.join(root, file_name)
                    open_yaml_file(file_path, csv_writer, root)

                    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_path", help="Path to the folder to scan for YAML files.")
    parser.add_argument("output_file", help="Path to the output CSV file.")
    args = parser.parse_args()
    
    process_yaml_files(args.folder_path, args.output_file)