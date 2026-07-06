import logging
import numpy as np
import pandas as pd
from itertools import product
from treelib import Tree

from BacTraq.helper import get_name, get_parent, increment_tuple_cluster
from BacTraq.reporting import _log_naming_events, _build_change_rows

logger = logging.getLogger('bactraq')

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
            elif get_parent(db_node_data) != parent:
                # Parent was renamed at a coarser threshold (e.g. merge/bridged cascade).
                # Fixing the stale history name whose prefix no longer matches.
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


def assign_cluster_database(data: pd.DataFrame, cluster_db: pd.DataFrame, change_log: list = None):
    ''' Get previous clonal complex cluster database.

    Parameters
    ----------
        - data (pd.DataFrame): The current cluster assignment from algorithm.
        - cluster_db (pd.DataFrame): The clonal complex cluster database. With index as sample id and columns as thresholds and value as cluster number.
        - change_log (list, optional): If a list is passed, per-level cluster change rows (old_cluster ->
          new_cluster edges — see `reporting._build_change_rows`) are appended to it in-place, for building
          the `--summary` table via `reporting.pivot_by_lineage`. Left untouched if None.

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

        history_samples = set(cluster_db.index)
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

                if change_log is not None:
                    change_log.extend(_build_change_rows(
                        overlap_cluster_, current_cluster_dict, db_cluster_dict, snp_d[t],
                        to_singleton, from_singleton, history_samples
                    ))

                #update them to new table
                replace_dict = overlap_cluster_[f'{snp_d[t]}_name_upd'].to_dict()
                df.loc[:, snp_d[t]:] = df.loc[:, snp_d[t]:].apply(upd_name_from_db, args=([replace_dict, t]), axis=1)
            else:
                logger.warning(f"No record of {snp_d[t]} cluster in history — skipping name tracking for this level")

        upd_db = df.sort_values(snp_d[0], na_position='last').copy()

        return df, upd_db
