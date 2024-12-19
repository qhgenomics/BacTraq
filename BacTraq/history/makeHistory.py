import pandas as pd
from datetime import datetime
import argparse
import os

def read_cluster_table(path, threshold: list):
    """ Read Cluster table from user input path csv file.

    Parameters:
        - path (str): path to cluster table csv file
        - threshold (list): list of snp threshold. Please make sure the threshold input in the same order as the columns in the cluster table
    
    Result:
        - pd.DataFrame: A dataframe with indexes as sample id, columns as SNP thresholds and values as cluster name
    """
    thres_names = [f'{t} SNPs' for t in threshold]
    sorted_thres = sorted(threshold, reverse=True)
    sorted_name = [f'{t} SNPs' for t in sorted_thres]

    df = pd.read_csv(path, header=0, index_col=0)
    df.columns = thres_names

    df = df.fillna(0)
    df[sorted_name[0]] = df[sorted_name[0]].astype('Int16')
    df = df.astype(str)
    # print(df.head())
    return df[sorted_name]

def row_str_to_tuple(row):
    values = row.values.tolist()
    cluster_tuples = []
    for v in values:
        if not v.startswith("0"):
            cluster_id = tuple(v.split('.'))
            cluster_id = tuple(map(int, cluster_id))
        else:
            cluster_id = None
        cluster_tuples.append(cluster_id)

    result = pd.Series(cluster_tuples, name = row.name, index=row.index)
    return result
            
def history_generate(path: str, thres_list: list, filename: str) -> None:
    df = read_cluster_table(path, thres_list)
    history = df.apply(row_str_to_tuple, axis = 1)
    # print(history.head())
    history.to_parquet(f'{filename}.parquet.gz', compression='gzip')
    print(f"History file: {filename}.parquet.gz")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Enter the path to your old cluster table. The file needs to be in comma separated format")
    parser.add_argument("-o", "--output", help="Save name file. Current support csv save file. If provided: `output.parquet.gz`, else: `inputfilename_history.parquet.gz")
    parser.add_argument("-t", "--threshold",  help="Enter your list of threshold. Make sure this match the column order in your cluster table. E.g: --threshold 20,10,5", required=True)
    
    args = parser.parse_args()
    threshold = args.threshold.split(',')
    threshold = list(map(int, threshold))
    his_file = args.file
    if args.output:
        output_name = args.output
    else:
        output_name = f"{os.path.basename(his_file)}_history"

    history_generate(his_file, thres_list=threshold, filename=output_name)

if __name__ == "__main__":
    main()

