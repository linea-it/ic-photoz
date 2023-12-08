import os
import pandas as pd
from sys import argv

def join_csv_files_from_directory(input_directory, output_file):
    csv_files_to_join = [os.path.join(input_directory, file) for file in os.listdir(input_directory) if file.endswith(".csv")]

    combined_data = pd.DataFrame()

    for file in csv_files_to_join:
        df = pd.read_csv(file)
        combined_data = pd.concat([combined_data, df], ignore_index=True)

    combined_data.to_csv(output_file, index=False)
    print(f"arquivo salvo {output_file}")

if __name__ == "__main__":
    output_csv_file = "sample_dp02.csv"
    input_directory = argv[1]

    join_csv_files_from_directory(input_directory, output_csv_file)