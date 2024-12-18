#!/usr/bin/env python
import os, sys
import pandas as pd
import numpy as np
import argparse
from sklearn.cluster import AgglomerativeClustering
import json
from itertools import product
from treelib import Tree, Node
from datetime import datetime
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

_DEFAULT_THRESHOLD = [20, 10, 5]

today = datetime.today().strftime("%Y%m%d")

# Tree supporter
def by_suffix(t):
    return t[:-1]
def by_last(t):
    return t[-1]
def get_name(node: tuple):
    return '.'.join(list(map(str, node)))
def get_parent(node):
    parent = node[:-1]
    if parent:
        return get_name(parent)
    else:
        return 'root'

def cluster_number_generate(df: pd.DataFrame, threshold: list) -> pd.DataFrame:
    ''' Reassign cluster number. Any cluster with only one element, their cluster number will be removed.
    
    Parameters 
    -----------
        - df (pd.DataFrame): A table with index as sample id and columns as assigned cluster numbers for different threshold.

        - threshold (list): Threshold to define a cluster. Can be a number or a list of threshold. Cannot be `NONE`.
            Default: [20, 10, 5]
            E.g: [10] or [20, 10, 5]

    Returns
    -------
        pd.DataFrame: A table with index as sample id and columns as assigned cluster numbers for different thresholds.

    '''
    def _genarate_name(row, snp_list) -> str:
        cluster = row[snp_list]
        clean = pd.Series(cluster).dropna().to_list()
        if len(clean) == len(snp_list):
            name = tuple(clean)
        else:
            name = None
        return name 
    
    def _new_reassign(row, groups, rename_dict):
        key = (row[groups].tolist())
        if len(key) == 1:
            if key[0] in rename_dict.keys():
                return rename_dict[key[0]]
            else:
                return None
        else:
            if tuple(key) in rename_dict.keys():
                return rename_dict[tuple(key)]
            else:
                None
    
    if threshold is None:
        raise AttributeError("Please provide a least on threshold in a list!")
        
    data = df.copy() # copy df

    rename = {}
    cols = [f'cluster{t}' for t in threshold]
    snp_d = [f'{t} SNPs' for t in threshold]
    for t in range(len(cols)):
        groups = cols[:t+1]
        df_groupby = pd.DataFrame(data.reset_index().groupby(groups)['sample'].size())
        df_groupby = df_groupby.replace(1, None)
        df_groupby = df_groupby.dropna(how='all')
        if t == 0:
            df_groupby.loc[:, 'sample'] = list(range(1, df_groupby.shape[0]+1))
            reassign = df_groupby['sample'].to_dict()
        else: 
            # df = df_groupby.reset_index(level=1)
            groups.append('reindex')
            df = df_groupby.set_index(df_groupby.groupby(level=t-1).cumcount().add(1), append=True)
            df.index.names = groups
            df = df.reset_index(level=t+1)
            reassign = df['reindex'].to_dict()

        data[snp_d[t]] = data.apply(_new_reassign, args=([cols[:t+1], reassign]), axis=1)
        data[snp_d[t]] = data[snp_d[t]].astype('Int8')
        naming = data.apply(_genarate_name, args=([snp_d[:t+1]]), axis=1)   
        rename[snp_d[t]] = naming
    cluster_name_df = pd.DataFrame(rename)    
    return cluster_name_df

def cluster_from_distance_matrix(path: str, threshold: list) -> pd.DataFrame:
    ''' Read in input distance matrix and generate clustering based on threshold
    
    Parameters 
    -----------
        - path (str): Path to distance matrix file. Tab separated. 
            Symetrix matix, have to have `Refernce` sample. 
            E.g:
                |         | Reference | Sample1 | Sample2 |
                |---------|-----------|---------|---------|
                |Reference|  0        |   20    |    12   |
                |Sample1  | 20        |    0    |    10   |

        - threshold (list): Threshold to define a cluster. Can be a number or a list of threshold. Cannot be `NONE`.
            Default: [20, 10, 5]
            E.g: [10] or [20, 10, 5]

    Returns
    -------
        - pd.DataFrame: A table with index as sample id and columns as assigned cluster numbers for different thresholds.
            E.g: 
            |         | 20 SNPs | 10 SNPs |  5 SNPs |
            |---------|---------|---------|---------|
            |Sample1  | (1, )   | (1, 1)  |(1, 1, 1)|
            |Sample2  | (2, )   | (2, 1)  |   None  |

    '''

    threshold = sorted(threshold, reverse=True)

    dist_mx = pd.read_csv(path, sep='\t', header=0, index_col=0)
    dist_mx.index.name = 'sample'
    if 'Reference' not in dist_mx.index:
        print("\033[91m {}\033[00m".format("Cannot find `Reference`. Make sure your SNPs dist matrix contains `Reference`"))
        sys.exit(1)

    clusters_assigned = {}
    for thres in threshold:
        clustering = AgglomerativeClustering(n_clusters=None, metric="precomputed", linkage="single", distance_threshold=thres+0.1) # +0.1 to include the threshold as cutoff
        clustering.fit(dist_mx)
        labels = clustering.labels_
        clusters_assigned[f'cluster{thres}'] = labels

    cluster_df = pd.DataFrame(clusters_assigned, index=dist_mx.index)
    cluster_df = cluster_df + 1 # add 1 so cluster number starting at 1
    cluster_df = cluster_df.drop(index='Reference')
    
    # re-assign number so that cluster with only one same won't be assign with any number
    cluster_number_reassign = cluster_number_generate(cluster_df, threshold)

    return cluster_number_reassign


def gather_cluster_database(db_path: str) -> pd.DataFrame:
    ''' Get previous clonal complex cluster database. 

    Parameters
    ----------
        - clonal_complex (str): Input the clonal complex to query cluster database.
    
    Return
    ------
        - pd.DataFrame: A table with index as sample ID and columns as cluster number at each threshold.
    '''
    cluster_db = json.loads(db_path)
    cluster_db = pd.read_json(cluster_db, dtype=False)

    return cluster_db

def new_in_first_list(list1: list, list2: list):
    ''' Return differences between 2 lists.

    Parameters
    -----------
        - list1 (list): A list of items.
        - list2 (list): A list of items.
    
    Return
    -------
        - different (list): a list of item in list1 that do not exist in list2.
    '''

    diff = list(set(list1) - set(list2))
    return diff

def intesection_betweet_cluster_dicts(query: dict, db: dict) -> pd.DataFrame:
    ''' Generate a dataframe showing the overlaps between clusters.

    Parameters
    ----------
        - query, db (dict): 2 dictionaries containing clusters information from this run and the db

    Return
    ------
        - pd.DataFrame: A dataframe with indexes as query clusters, columns as db cluster and values as if they are overlap (True) or not (False)
    '''
    intersection = pd.DataFrame(index = query.keys(), columns = db.keys(), dtype="boolean")
    for right, left in product(query, db):
        intersection.loc[right, left] = bool(list(set(query[right]) &  set(db[left])))
    return intersection

def match_db_cluster(row, db_cluster):
        mask = row.values
        clusters = np.empty(len(db_cluster), dtype=object)
        clusters[:] = db_cluster
        overlap = clusters[mask]
        if len(overlap) > 0:
            return overlap
        else:
            return None
        
def increment_tuple_cluster(cluster_numb: tuple, add_on: int) -> list:
        name = list(cluster_numb)
        name[-1] = cluster_numb[-1] + 1
        return tuple(name)          

def assign_cluster_from_db(row, snp_t, old_t):
    overlap = row.loc[snp_t].item()
    node = row.name
    parent = get_parent(node)
    db_node_siblings = old_t.is_branch(parent)
    if db_node_siblings:
        siblings_info = {x :old_t.get_node(x).data for x in db_node_siblings}
        big_sis_info = sorted(siblings_info.values(), reverse=True)[0]
    else: 
        # print('this is entire new case')
        new_cluster = row.name
        old_t.create_node(get_name(new_cluster), get_name(new_cluster), parent=parent, data=new_cluster)
        return new_cluster
    if overlap is not None and len(list(overlap)) == 1:
        # print('Only overlap 1 cluster in db which is', overlap.item())
        db_node_data = list(overlap)[0]
        # print(row['split'], type(row.split))
        # print(np.all(row['split']))
        if np.all(row['split']): 
            # print('This is split case')
            new_cluster = increment_tuple_cluster(big_sis_info, 1)
            old_t.create_node(get_name(new_cluster), get_name(new_cluster), parent=parent, data=new_cluster)
            # cluster_name = get_name(new_cluster)
        else:
            # print('This is expand case')
            new_cluster = db_node_data
    else:

        new_cluster = increment_tuple_cluster(big_sis_info, 1)
        old_t.create_node(get_name(new_cluster), get_name(new_cluster), parent=parent, data=new_cluster)
    
    return new_cluster

def upd_name_from_db(row, name_change, index):
    new_values = []
    for cluster in row.values.tolist():
        if cluster is None:
            new_values.append(cluster)
        elif cluster[:index+1] != name_change[cluster[:index+1]]:
            rename = list(name_change[cluster[:index+1]])
            cluster = list(cluster)
            cluster[:index+1] = rename
            cluster = tuple(cluster)
            new_values.append(cluster)
        else:
            new_values.append(cluster)
    pd_value = pd.Series(new_values, index=row.index, name=row.name)
    return pd_value

def assign_cluster_database(data: pd.DataFrame, cluster_db: pd.DataFrame):
    ''' Get previous clonal complex cluster database. 

    Parameters
    ----------
        - cluster_db (pd.DataFrame): The clonal complex cluster database. With index as sample id and columns as thresholds and value as cluster number.
        - cluster_db (pd.DataFrame): The clonal complex cluster database. With index as sample id and columns as thresholds and value as cluster number.
    
    Return
    ------
        - pd.DataFrame: A table with index as sample ID and columns as cluster number reassigned based on cluster database.

    TODO:
        - Update database if new thresholds were used 
        - Assign number or Update cluster structure based on previous run (search for db and date created)
        - Return in a table.
    '''

    if cluster_db is None:
        cluster_db = data.copy()
        return data, cluster_db

    # check new threshold
    new_threshold = new_in_first_list(data.columns, cluster_db.columns)
    
    # updata new threshold to database 
    if new_threshold: 
        upd_db = data.copy()
        upd_db = upd_db.dropna(how='all')
        return data, upd_db
    else:
        old_tree = Tree()
        old_tree.create_node('Root', 'root')

        df = data.copy()
        snp_d = list(data.columns)
        for t in range(len(snp_d)):
            db_cluster_dict = cluster_db.reset_index().groupby(snp_d[t])['sample'].agg(list).to_dict()
            current_cluster_dict = df.reset_index().groupby(snp_d[t])['sample'].agg(list).to_dict()

            db_cluster = list(db_cluster_dict.keys())

            if db_cluster:
            # current_cluster = list(cc152_sep_dict.keys())

                for n in db_cluster:
                    old_tree.create_node(get_name(n), get_name(n), parent=get_parent(n), data = n)
                
                #create overlap between new and old cluster then reassign
                
                
                overlap_cluster_ = intesection_betweet_cluster_dicts(current_cluster_dict, db_cluster_dict) 
                overlap_cluster_[f'{snp_d[t]}_db_match'] = overlap_cluster_.apply(match_db_cluster, args=([list(db_cluster_dict.keys())]), axis=1)
                overlap_cluster_['split'] = overlap_cluster_.loc[~overlap_cluster_[(f'{snp_d[t]}_db_match',)].isna(), (f'{snp_d[t]}_db_match',)].duplicated(keep=False)
                overlap_cluster_['split'] = overlap_cluster_['split'].fillna(False)
                overlap_cluster_[f'{snp_d[t]}_name_upd'] = overlap_cluster_.apply(assign_cluster_from_db, args=([f'{snp_d[t]}_db_match', old_tree]), axis=1)
                
                #update them to new table
                replace_dict = overlap_cluster_.loc[:, (f'{snp_d[t]}_name_upd',)].to_dict()
                df.loc[:, snp_d[t]:] = df.loc[:, snp_d[t]:].apply(upd_name_from_db, args=([replace_dict, t]), axis=1)
                overlap_cluster_.to_csv(f"{snp_d[t]}_cc152_pipeline.csv")
            else:
                print(f"no record of {snp_d[t]} cluster in history_input")
    
        upd_db = df.sort_values(snp_d[0]).dropna(how='all').copy()

        return df, upd_db

def pretty_name_save(data: pd.DataFrame, matrix_file: str, name: str) -> None: 
    df = pd.DataFrame(index = data.index)
    input_file = os.path.basename(matrix_file)
    for c in data.columns:
        df[c] = data[c].dropna(how='all').apply(lambda x: '.'.join(list(map(str, x))))
    df.to_csv(name)
    print(f'Clusters for {input_file} is saved in {name}')
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snpdist_matrix", help="Enter the path to SNPdist matrix.")
    parser.add_argument("output", help="Save name file. Current support csv save file")
    parser.add_argument("-t", "--threshold",  help="Enter your list of threshold. E.g: --threshold 20,10,5. Default: 20,10,5")
    parser.add_argument("--history",  help="Enter the path to previous/databases clusters to match name. This should include sample id and clusters at chosen thresholds. If empty then will have this run as history.")
    args = parser.parse_args()

    distance_matrix = args.snpdist_matrix
    save_name = args.output
    if args.threshold:
        snp_thres = args.threshold.split(',')
        snp_thres = list(map(int, snp_thres))
        snp_thres = sorted(snp_thres, reverse=True)
    else:
        snp_thres = _DEFAULT_THRESHOLD
    
    if args.history:
        db_cluster = pd.read_parquet(args.history)
        for c in db_cluster.columns:
            db_cluster.loc[:, c] = db_cluster.loc[:, c].apply(lambda x: tuple(x) if x is not None else None)
    else: 
        db_cluster = None
    
    this_run_cluster = cluster_from_distance_matrix(distance_matrix, snp_thres)

    rename_with_db, new_db = assign_cluster_database(this_run_cluster, cluster_db=db_cluster)

    pretty_name_save(rename_with_db, distance_matrix, args.output)
    new_db.to_parquet(f'{today}_cluster_history.parquet.gz', compression='gzip')
    print(f"New history cluster is created for this run and saved in `{today}_cluster.parquet.gz` ")

if __name__=="__main__":
    main()

    


    

