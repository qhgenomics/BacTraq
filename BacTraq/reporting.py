import os
import logging
import pandas as pd

from BacTraq.helper import get_name, get_parent

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
                former_singletons = sorted(s for s in members if s in from_singleton)
                if get_parent(new_cluster) != get_parent(db_cluster):
                    # Parent cluster was renamed at a coarser threshold (merge/bridged cascade).
                    logger.info(
                        f"[{snp_col}] CASCADED RENAME  cluster {get_name(db_cluster)} "
                        f"-> {new_name} (parent renamed at coarser threshold); "
                        f"members: [{', '.join(members)}]"
                    )
                else:
                    # Consistent match but singletons joined, triggering a rename (bridged reclassification)
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


def _change_row(snp_col: str, old_cluster: str, new_cluster: str, event: str,
                 samples: list, n_old=None, n_new=None, reason: str = '') -> dict:
    ''' Build one edge-list row for the cluster change table: an old_cluster -> new_cluster
    flow of `samples`, tagged with the atomic event that produced it. One (old, new) pair
    can appear multiple times across a run (e.g. a merge produces one row per contributing
    history cluster) — this is intentionally an edge list, not one row per cluster. '''
    samples = sorted(samples)
    return {
        'level': snp_col,
        'old_cluster': old_cluster,
        'new_cluster': new_cluster,
        'event': event,
        'reason': reason,
        'n_old': n_old,
        'n_new': n_new,
        'n_shared': len(samples),
        'samples': ', '.join(samples),
    }


def _build_change_rows(overlap_df: pd.DataFrame, current_cluster_dict: dict,
                        db_cluster_dict: dict, snp_col: str,
                        to_singleton: dict, from_singleton: dict,
                        history_samples: set) -> list:
    ''' Turn one threshold level's naming-event data into cluster-to-cluster edge rows
    (old_cluster -> new_cluster, weighted by shared sample count), for the --summary table.
    'SINGLETON' and 'NEW-SAMPLE' are pseudo-nodes standing in for "unclustered" and "sample
    absent from history", respectively.

    For CONTINUED edges (single, non-split match) whose name changed, 'reason' explains why:
    a rename with no membership change at this level is either inherited from a coarser-level
    rename of the parent cluster (cascade), or caused by a singleton joining at this level.

    Parameters
    ----------
        - overlap_df: per-current-cluster overlap table for this level (from assign_cluster_database)
        - current_cluster_dict, db_cluster_dict: cluster tuple -> member sample list (current / history)
        - snp_col: threshold column being examined (e.g. '20 SNPs')
        - to_singleton, from_singleton: sample -> cluster tuple, from detect_singleton_transitions
        - history_samples: every sample id present anywhere in the history db

    Returns
    -------
        - list[dict]: rows for the cluster change table, see `_change_row`
    '''
    name_upd_col = f'{snp_col}_name_upd'
    db_match_col = f'{snp_col}_db_match'
    rows = []

    for curr_cluster, row in overlap_df.iterrows():
        new_cluster = row[name_upd_col]
        match = row[db_match_col]
        is_split = bool(row['split'])
        members = set(current_cluster_dict.get(curr_cluster, []))
        new_name = get_name(new_cluster)
        n_new = len(members)
        accounted = set()

        if match is not None:
            for db_cluster in match:
                accounted |= (members & set(db_cluster_dict.get(db_cluster, [])))

        leftover = members - accounted
        joined_from_singleton = {m for m in leftover if m in from_singleton}
        brand_new = {m for m in leftover - joined_from_singleton if m not in history_samples}
        unexplained = leftover - joined_from_singleton - brand_new

        if match is not None:
            is_merge = len(list(match)) > 1
            event = 'MERGE' if is_merge else ('SPLIT' if is_split else 'CONTINUED')
            for db_cluster in match:
                db_members = set(db_cluster_dict.get(db_cluster, []))
                shared = members & db_members
                if not shared:
                    continue
                reason = ''
                if event == 'CONTINUED' and new_name != get_name(db_cluster):
                    if joined_from_singleton:
                        reason = 'bridged: singleton joined at this level'
                    elif get_parent(new_cluster) != get_parent(db_cluster):
                        reason = 'cascaded: parent cluster renamed at a coarser threshold'
                    else:
                        reason = 'renamed (unexplained)'
                rows.append(_change_row(snp_col, get_name(db_cluster), new_name, event,
                                         shared, n_old=len(db_members), n_new=n_new, reason=reason))

        if joined_from_singleton:
            rows.append(_change_row(snp_col, 'SINGLETON', new_name, 'FROM-SINGLETON',
                                     joined_from_singleton, n_new=n_new))
        if brand_new:
            rows.append(_change_row(snp_col, 'NEW-SAMPLE', new_name, 'NEW',
                                     brand_new, n_new=n_new))
        if unexplained:
            # Shouldn't normally happen; kept so no sample silently disappears from the table.
            rows.append(_change_row(snp_col, 'UNKNOWN', new_name, 'NEW', unexplained, n_new=n_new))

    # Departures to singleton, grouped by the history cluster they left.
    departures = {}
    for sample, prev_cluster in to_singleton.items():
        departures.setdefault(prev_cluster, []).append(sample)
    for prev_cluster, samples in departures.items():
        rows.append(_change_row(snp_col, get_name(prev_cluster), 'SINGLETON', 'TO-SINGLETON',
                                 samples, n_old=len(db_cluster_dict.get(prev_cluster, []))))

    return rows


def pivot_by_lineage(changes_df: pd.DataFrame, level_order: list,
                    current_assignments: pd.DataFrame = None,
                    history_assignments: pd.DataFrame = None) -> pd.DataFrame:
    ''' Build the --summary cluster-change table: one row per history cluster lineage, with a
    Status/New Name column pair per SNP threshold, plus a Members column for the finest level.

    Parameters
    ----------
        - changes_df (pd.DataFrame): Edge-list change table from `_build_change_rows`.
        - level_order (list): SNP threshold column names, coarsest first.
        - current_assignments (pd.DataFrame, optional): This run's final per-sample cluster table.
          If given, singleton samples at the finest level get their own trailing rows.
        - history_assignments (pd.DataFrame, optional): Same-shaped history cluster table, used to
          tell an unchanged singleton apart from one that just left a cluster.

    Return
    ------
        - pd.DataFrame: Indexed by lineage (or sample id for singleton rows), one Status/New Name
          pair per level. Status is Unchanged / Renamed / Expanded / Reduced / Split / Merged /
          Dissolved / "-" for lineage rows, or Clustered / Singleton / Departed for singleton rows.
    '''
    df = changes_df.copy()
    finest_level = level_order[-1]

    def lineage_of(name):
        if name in ('SINGLETON', 'NEW-SAMPLE', 'UNKNOWN'):
            return None
        return name.split('.')[0]

    df['lineage'] = df['old_cluster'].apply(lineage_of)
    lineages = sorted(df.loc[df['lineage'].notna(), 'lineage'].unique(), key=int)

    rows = []
    for lineage in lineages:
        row = {'old_cluster': lineage}
        finest_targets = []
        for level in level_order:
            sub = df[(df['lineage'] == lineage) & (df['level'] == level)]
            outbound = sub[sub['event'].isin(['CONTINUED', 'MERGE', 'SPLIT'])]
            dissolved = sub[sub['event'] == 'TO-SINGLETON']
            new_names = sorted(outbound['new_cluster'].unique())

            if sub.empty:
                status, new_name = '-', '-'
            elif not new_names and not dissolved.empty:
                status, new_name = 'Dissolved', '-'
            elif len(new_names) > 1:
                status = 'Split' + (' (+singleton departures)' if not dissolved.empty else '')
                new_name = '/'.join(new_names)
            else:
                target = new_names[0]
                is_merge = (outbound['event'] == 'MERGE').any()
                inflow = df[(df['level'] == level) & (df['new_cluster'] == target) &
                            (df['event'].isin(['FROM-SINGLETON', 'NEW']))]
                old_name = outbound.iloc[0]['old_cluster']
                if is_merge:
                    status = 'Merged'
                elif not dissolved.empty:
                    status = 'Reduced'
                elif not inflow.empty:
                    status = 'Expanded'
                elif target != old_name:
                    status = 'Renamed'
                else:
                    status = 'Unchanged'
                new_name = target

            if level == finest_level:
                finest_targets = new_names

            row[(level, 'Status')] = status
            row[(level, 'New Name')] = new_name

        # One member group per finest-level target, in the same order as New Name, so a
        # split's "12.1.1/9.1.1" lines up with "s8, s9 / n1, s10" group-for-group.
        member_groups = []
        for target in finest_targets:
            hits = df[(df['level'] == finest_level) & (df['new_cluster'] == target)]
            target_members = set()
            for s in hits['samples']:
                if s:
                    target_members.update(s.split(', '))
            member_groups.append(', '.join(sorted(target_members)) if target_members else '-')
        row[('Members', f'({finest_level})')] = ' / '.join(member_groups) if member_groups else '-'

        rows.append(row)
    # tracing singleton events
    if current_assignments is not None:
        is_singleton = current_assignments[finest_level].apply(lambda v: v is None)
        finest_singletons = current_assignments.index[is_singleton]
        for sample in sorted(finest_singletons):
            srow = {'old_cluster': sample}
            for level in level_order:
                val = current_assignments.loc[sample, level]
                if val is not None:
                    srow[(level, 'Status')] = 'Clustered'
                    srow[(level, 'New Name')] = get_name(val)
                    continue

                prev = None
                if history_assignments is not None and sample in history_assignments.index:
                    prev = history_assignments.loc[sample, level]

                if prev is not None:
                    srow[(level, 'Status')] = f'Departed left {get_name(prev)}'
                    srow[(level, 'New Name')] = '-'
                else:
                    srow[(level, 'Status')] = 'Singleton'
                    srow[(level, 'New Name')] = '-'
            srow[('Members', f'({finest_level})')] = sample
            rows.append(srow)

    pivot = pd.DataFrame(rows).set_index('old_cluster')
    pivot.columns = pd.MultiIndex.from_tuples(pivot.columns)
    pivot.index.name = 'Old'
    return pivot


def pretty_name_save(data: pd.DataFrame, matrix_file: str, name: str) -> None:
    df = pd.DataFrame(index = data.index)
    input_file = os.path.basename(matrix_file)
    for c in data.columns:
        df[c] = data[c].dropna(how='all').apply(lambda x: '.'.join(list(map(str, x))))
    df.to_csv(name)
    logger.info(f"Clusters for {input_file} saved to {name}")
    return None
