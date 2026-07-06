import sys
import logging
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

logger = logging.getLogger('bactraq')

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

    data = df.copy()

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
