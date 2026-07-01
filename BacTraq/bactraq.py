#!/usr/bin/env python
import os, sys
import logging
import pandas as pd
import numpy as np
import argparse
from sklearn.cluster import AgglomerativeClustering
from itertools import product
from treelib import Tree
from datetime import datetime
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

logger = logging.getLogger('bactraq')

def _setup_logging(log_path: str) -> None:
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(fmt)

    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    logger.addHandler(console)
    logger.addHandler(fh)

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

def cluster_from_distance_matrix(path: str, threshold: list, ref_exist: bool) -> pd.DataFrame:
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
    if 'Reference' not in dist_mx.index and ref_exist:
        logger.error("Cannot find `Reference`. Make sure your SNPs dist matrix contains `Reference`")
        sys.exit(1)

    clusters_assigned = {}
    for thres in threshold:
        clustering = AgglomerativeClustering(n_clusters=None, metric="precomputed", linkage="single", distance_threshold=thres+0.1) # +0.1 to include the threshold as cutoff
        clustering.fit(dist_mx)
        labels = clustering.labels_
        clusters_assigned[f'cluster{thres}'] = labels

    cluster_df = pd.DataFrame(clusters_assigned, index=dist_mx.index)
    cluster_df = cluster_df + 1 # add 1 so cluster number starting at 1
    if ref_exist:
        cluster_df = cluster_df.drop(index='Reference')
    
    # re-assign number so that cluster with only one same won't be assign with any number
    cluster_number_reassign = cluster_number_generate(cluster_df, threshold)

    return cluster_number_reassign

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
        
def increment_tuple_cluster(cluster_numb: tuple) -> list:
        name = list(cluster_numb)
        name[-1] = cluster_numb[-1] + 1
        return tuple(name)


def detect_singleton_transitions(df: pd.DataFrame, cluster_db: pd.DataFrame, snp_col: str):
    ''' Identify samples that moved to or from singleton status between runs.

    Parameters
    ----------
        df         : current-run cluster assignments (None = singleton at this level)
        cluster_db : history cluster assignments
        snp_col    : threshold column being examined (e.g. '20 SNPs')

    Returns
    -------
        to_singleton             : {sample: prev_cluster_tuple}  — was named, now singleton
        from_singleton           : {sample: curr_cluster_tuple}  — was singleton, now named
        affected_history_clusters: set of history cluster tuples that lost ≥1 member to singleton
    '''
    shared = df.index.intersection(cluster_db.index)
    to_singleton = {}
    from_singleton = {}

    for sample in shared:
        prev = cluster_db.loc[sample, snp_col]
        curr = df.loc[sample, snp_col]
        if prev is not None and curr is None:
            to_singleton[sample] = prev
        elif prev is None and curr is not None:
            from_singleton[sample] = curr

    return to_singleton, from_singleton, set(to_singleton.values())


def _inject_split_if_affected(row, db_match_col: str, affected_clusters: set) -> bool:
    ''' Return True if this current cluster should be treated as a split because its
    matching history cluster lost members to singleton status this run. '''
    if np.all(row['split']):
        return True
    match = row[db_match_col]
    return match is not None and any(m in affected_clusters for m in match)


def _log_naming_events(overlap_df: pd.DataFrame, current_cluster_dict: dict,
                        db_cluster_dict: dict, snp_col: str,
                        to_singleton: dict, from_singleton: dict,
                        affected_history_clusters: set) -> None:
    ''' Log every cluster naming event for one threshold level after names are assigned. '''
    name_upd_col = f'{snp_col}_name_upd'
    db_match_col = f'{snp_col}_db_match'
    logged_affected = set()

    for curr_cluster, row in overlap_df.iterrows():
        new_cluster = row[name_upd_col]
        match = row[db_match_col]
        is_split = bool(row['split'])
        members = sorted(current_cluster_dict.get(curr_cluster, []))
        new_name = get_name(new_cluster)

        if match is None:
            # No overlap with any history cluster
            former_singletons = sorted(s for s in members if s in from_singleton)
            if former_singletons:
                logger.info(
                    f"[{snp_col}] NEW (from singletons)  cluster {new_name}: "
                    f"previously-singleton [{', '.join(former_singletons)}] formed a new cluster; "
                    f"all members: [{', '.join(members)}]"
                )
            else:
                logger.info(
                    f"[{snp_col}] NEW  cluster {new_name}: [{', '.join(members)}]"
                )

        elif len(list(match)) > 1:
            # Multiple history clusters merged into this one
            merged_names = ', '.join(get_name(c) for c in match)
            logger.info(
                f"[{snp_col}] MERGE  clusters [{merged_names}] -> "
                f"new cluster {new_name}; members: [{', '.join(members)}]"
            )

        else:
            db_cluster = list(match)[0]
            db_members = sorted(db_cluster_dict.get(db_cluster, []))

            if db_cluster in affected_history_clusters:
                # Gap 2: some members of this history cluster became singletons this run
                departed = sorted(s for s, c in to_singleton.items() if c == db_cluster)
                logged_affected.add(db_cluster)
                logger.info(
                    f"[{snp_col}] SPLIT-TO-SINGLETON  "
                    f"[{', '.join(departed)}] left cluster {get_name(db_cluster)} and became singleton(s); "
                    f"remaining members [{', '.join(members)}] renamed: "
                    f"{get_name(db_cluster)} -> {new_name}"
                )

            elif is_split:
                logger.info(
                    f"[{snp_col}] SPLIT  cluster {get_name(db_cluster)} "
                    f"(had [{', '.join(db_members)}]) -> "
                    f"new cluster {new_name}: [{', '.join(members)}]"
                )

            elif new_cluster != db_cluster:
                # Consistent match but singletons joined, triggering a rename (bridged reclassification)
                former_singletons = sorted(s for s in members if s in from_singleton)
                logger.info(
                    f"[{snp_col}] BRIDGED RECLASSIFICATION  "
                    f"[{', '.join(former_singletons)}] (previously singleton) "
                    f"joined cluster {get_name(db_cluster)}; "
                    f"renamed: {get_name(db_cluster)} -> {new_name}; "
                    f"all members: [{', '.join(members)}]"
                )

            else:
                # Name retained
                new_members = sorted(s for s in members if s not in db_members)
                msg = (
                    f"[{snp_col}] CONSISTENT  cluster {get_name(db_cluster)}: "
                    f"name retained as {new_name}; members: [{', '.join(members)}]"
                )
                if new_members:
                    msg += f"; new samples added: [{', '.join(new_members)}]"
                logger.info(msg)

    # History clusters where every member became a singleton — no current cluster maps to them
    for hist_cluster in affected_history_clusters - logged_affected:
        departed = sorted(s for s, c in to_singleton.items() if c == hist_cluster)
        logger.info(
            f"[{snp_col}] DISSOLVED  cluster {get_name(hist_cluster)}: "
            f"all members [{', '.join(departed)}] became singletons"
        )


def assign_cluster_from_db(row, snp_t, old_t, current_cluster_dict=None, singleton_samples=None):
    overlap = row[snp_t]
    node = row.name
    parent = get_parent(node)
    db_node_siblings = old_t.is_branch(parent)
    if db_node_siblings:
        siblings_info = {x: old_t.get_node(x).data for x in db_node_siblings}
        big_sis_info = sorted(siblings_info.values(), reverse=True)[0]
    else:
        new_cluster = row.name
        old_t.create_node(get_name(new_cluster), get_name(new_cluster), parent=parent, data=new_cluster)
        return new_cluster
    if overlap is not None and len(list(overlap)) == 1:
        db_node_data = list(overlap)[0]
        if np.all(row['split']):
            new_cluster = increment_tuple_cluster(big_sis_info)
            old_t.create_node(get_name(new_cluster), get_name(new_cluster), parent=parent, data=new_cluster)
        else:
            # Check for bridged reclassification: previously-singleton samples now in this cluster
            bridged = set()
            if current_cluster_dict is not None and singleton_samples is not None:
                current_samples = set(current_cluster_dict.get(node, []))
                bridged = current_samples & singleton_samples
            if bridged:
                new_cluster = increment_tuple_cluster(big_sis_info)
                old_t.create_node(get_name(new_cluster), get_name(new_cluster), parent=parent, data=new_cluster)
            else:
                new_cluster = db_node_data
    else:
        new_cluster = increment_tuple_cluster(big_sis_info)
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
        - data (pd.DataFrame): The current cluster assignment from algorithm.
        - cluster_db (pd.DataFrame): The clonal complex cluster database. With index as sample id and columns as thresholds and value as cluster number.
    
    Return
    ------
        - pd.DataFrame: A table with index as sample ID and columns as cluster number reassigned based on cluster database.
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
                for n in db_cluster:
                    old_tree.create_node(get_name(n), get_name(n), parent=get_parent(n), data=n)

                # Detect singleton transitions before pairwise comparison
                to_singleton, from_singleton, affected_history_clusters = detect_singleton_transitions(
                    df, cluster_db, snp_d[t]
                )

                singleton_samples = set(cluster_db.index[cluster_db[snp_d[t]].isna()])
                overlap_cluster_ = intesection_betweet_cluster_dicts(current_cluster_dict, db_cluster_dict)
                # Rename tuple column labels to strings so pandas doesn't create a MultiIndex,
                # which would cause row[string_col] to return a Series instead of a scalar.
                overlap_cluster_.columns = [get_name(c) for c in overlap_cluster_.columns]
                overlap_cluster_[f'{snp_d[t]}_db_match'] = overlap_cluster_.apply(match_db_cluster, args=([list(db_cluster_dict.keys())]), axis=1)
                overlap_cluster_['split'] = overlap_cluster_.loc[~overlap_cluster_[f'{snp_d[t]}_db_match'].isna(), f'{snp_d[t]}_db_match'].duplicated(keep=False)
                overlap_cluster_['split'] = overlap_cluster_['split'].fillna(False)

                # Gap 2: treat current clusters as split if their history cluster lost members to singleton
                if affected_history_clusters:
                    overlap_cluster_['split'] = overlap_cluster_.apply(
                        _inject_split_if_affected,
                        args=([f'{snp_d[t]}_db_match', affected_history_clusters]),
                        axis=1
                    )

                overlap_cluster_[f'{snp_d[t]}_name_upd'] = overlap_cluster_.apply(
                    assign_cluster_from_db,
                    args=([f'{snp_d[t]}_db_match', old_tree, current_cluster_dict, singleton_samples]),
                    axis=1
                )

                _log_naming_events(
                    overlap_cluster_, current_cluster_dict, db_cluster_dict, snp_d[t],
                    to_singleton, from_singleton, affected_history_clusters
                )

                #update them to new table
                replace_dict = overlap_cluster_[f'{snp_d[t]}_name_upd'].to_dict()
                df.loc[:, snp_d[t]:] = df.loc[:, snp_d[t]:].apply(upd_name_from_db, args=([replace_dict, t]), axis=1)
            else:
                logger.warning(f"No record of {snp_d[t]} cluster in history — skipping name tracking for this level")
    
        upd_db = df.sort_values(snp_d[0], na_position='last').copy()

        return df, upd_db

def pretty_name_save(data: pd.DataFrame, matrix_file: str, name: str) -> None: 
    df = pd.DataFrame(index = data.index)
    input_file = os.path.basename(matrix_file)
    for c in data.columns:
        df[c] = data[c].dropna(how='all').apply(lambda x: '.'.join(list(map(str, x))))
    df.to_csv(name)
    logger.info(f"Clusters for {input_file} saved to {name}")
    return None

def main():
    parser = argparse.ArgumentParser(
        description="SNP-distance clustering with consistent naming across runs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-d", "--distance-matrix", required=True, metavar="FILE",
                        help="Path to SNPdist matrix (tab-separated).")
    parser.add_argument("-o", "--output", required=True, metavar="FILENAME",
                        help="Output base name. Produces <output>.csv and <output>_history.parquet.gz.")
    parser.add_argument("-t", "--threshold", metavar="THRESHOLD",
                        help="Comma-separated SNP thresholds, largest first.\nE.g: --threshold 20,10,5\nDefault: 20,10,5")
    parser.add_argument("--history", metavar="HISTORY",
                        help="Path to history .parquet.gz from a previous run.\nGenerate with `bactraq-history`. Omit for clustering only.")
    parser.add_argument("--nRef", action='store_false',
                        help="No Reference sample in the SNP distance matrix.")
    args = parser.parse_args()

    distance_matrix = args.distance_matrix
    save_name = f"{args.output}.csv"
    log_file = f"{args.output}_bactraq.log"
    _setup_logging(log_file)
    logger.info(f"BacTraq run started — input: {distance_matrix}")

    if args.threshold:
        snp_thres = sorted(map(int, args.threshold.split(',')), reverse=True)
    else:
        snp_thres = _DEFAULT_THRESHOLD
    logger.info(f"Thresholds: {snp_thres}")

    if args.history:
        logger.info(f"History file: {args.history}")
        db_cluster = pd.read_parquet(args.history)
        db_cluster = db_cluster.apply(
            lambda col: col.map(lambda x: tuple(x) if x is not None else None)
        )
    else:
        logger.info("No history file provided — clustering only, no name tracking")
        db_cluster = None

    this_run_cluster = cluster_from_distance_matrix(distance_matrix, snp_thres, args.nRef)
    logger.info(f"Clustering complete — {len(this_run_cluster)} samples")

    rename_with_db, new_db = assign_cluster_database(this_run_cluster, cluster_db=db_cluster)

    pretty_name_save(rename_with_db, distance_matrix, save_name)

    new_history_file = f'{args.output}_history.parquet.gz'
    new_db.to_parquet(new_history_file, compression='gzip')
    logger.info(f"History updated and saved to {new_history_file}")
    logger.info(f"Log saved to {log_file}")

if __name__=="__main__":
    main()