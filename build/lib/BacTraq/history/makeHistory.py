import pandas as pd
from datetime import datetime
import argparse

def read_cluster_table(path, threshold: list):
    """ Read Cluster table from user input path csv file.

    Parameters:
        - path (str): path to cluster table csv file
        - threshold (list): list of snp threshold. Please make sure the threshold input in the same order as the columns in the cluster table
    
    Result:
        - pd.DataFrame: A dataframe with indexes as sample id, columns as SNP thresholds and values as cluster name
    """
    thres_names = [f'{t} SMPs' for t in threshold]
    sorted_thres = sorted(threshold, reversed=True)
    sorted_name = [f'{t} SNPs' for t in sorted_thres ]

    df = pd.read_csv(path, head=0, index_col=0)
    df.columns = thres_names

    return df[[sorted_name]]

def row_str_to_tuple(row):
    values = row.values.tolist()
    cluster_tuples = []
    for v in values:
        if v:
            cluster_id = tuple(v.split('.'))
        else:
            cluster_id = None
        cluster_tuples.append(cluster_id)

    result = pd.Series(cluster_tuples, name = row.name, index=row.index)
    return result
            
def history_generate(path: str, thres_list: list, filename: str = None) -> None:
    df = read_cluster_table(path, thres_list)
    history = df.apply(row_str_to_tuple, axis = 1)
    if filename:
        history.to_parquet(f'{filename}.parquet.gz', compress='gzip')
    else:
        today = datetime.today().strftime("%Y%m%d")
        history.to_parquet(f'{today}_history_cluster.parquet.gz', compression='gzip')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Enter the path to SNPdist matrix.")
    parser.add_argument("-o", "--output", help="Save name file. Current support csv save file")
    parser.add_argument("-t", "--threshold",  help="Enter your list of threshold. E.g: --threshold 20,10,5.", required=True)
    
    args = parser.parse_args()
    threshold = args.threshold.split(',')
    his_file = args.file
    output_name = args.output

    history_generate(his_file, thres_list=threshold, filename=output_name)

