''' Small pure helpers for working with hierarchical cluster-name tuples, e.g. (1, 2, 1) '''

def by_suffix(t: tuple) -> tuple:
    return t[:-1]


def by_last(t: tuple):
    return t[-1]


def get_name(node: tuple) -> str:
    return '.'.join(list(map(str, node)))


def get_parent(node: tuple) -> str:
    parent = node[:-1]
    if parent:
        return get_name(parent)
    else:
        return 'root'


def increment_tuple_cluster(cluster_numb: tuple) -> tuple:
    name = list(cluster_numb)
    name[-1] = cluster_numb[-1] + 1
    return tuple(name)
