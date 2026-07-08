import logging
import argparse
import pandas as pd

from BacTraq.clustering import cluster_from_distance_matrix
from BacTraq.cluster_trace import assign_cluster_database
from BacTraq.reporting import _setup_logging, pretty_name_save, pivot_by_lineage, build_rename_trace

logger = logging.getLogger('bactraq')

_DEFAULT_THRESHOLD = [20, 10, 5]

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
                        default=','.join(map(str, _DEFAULT_THRESHOLD)),
                        help="Comma-separated SNP thresholds, largest first.\nE.g: --threshold 20,10,5")
    parser.add_argument("--history", metavar="HISTORY",
                        help="Path to history .parquet.gz from a previous run.\nGenerate with `bactraq-history`. Omit for clustering only.")
    parser.add_argument("--nRef", action='store_false',
                        help="No Reference sample in the SNP distance matrix. `Reference` is usually included when you run the workflow with Snippy and Snippy-core.")
    parser.add_argument("--summary", action='store_true',
                        help="Also write <output>_summary.xlsx: an 'Events' sheet (one row per history "
                             "cluster lineage, Status/New Name per threshold, plus singleton samples) and "
                             "a 'Rename Trace' sheet (one row per sample, Old/New cluster name per "
                             "threshold) comparing this run against --history. No-op without --history.")
    args = parser.parse_args()

    distance_matrix = args.distance_matrix
    save_name = f"{args.output}.csv"
    log_file = f"{args.output}_bactraq.log"
    _setup_logging(log_file)
    logger.info(f"BacTraq run started — input: {distance_matrix}")

    snp_thres = sorted(map(int, args.threshold.split(',')), reverse=True)
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

    change_log = [] if args.summary else None
    rename_with_db, new_db = assign_cluster_database(this_run_cluster, cluster_db=db_cluster, change_log=change_log)

    pretty_name_save(rename_with_db, distance_matrix, save_name)

    if args.summary:
        if change_log:
            changes_df = pd.DataFrame(change_log)
            level_order = list(this_run_cluster.columns)  # coarsest -> finest
            events_table = pivot_by_lineage(changes_df, level_order, current_assignments=rename_with_db, history_assignments=db_cluster)
            rename_trace_table = build_rename_trace(rename_with_db, db_cluster, level_order)
            summary_file = f"{args.output}_summary.xlsx"
            with pd.ExcelWriter(summary_file) as writer:
                events_table.to_excel(writer, sheet_name='Events')
                rename_trace_table.to_excel(writer, sheet_name='Rename Trace')
            logger.info(f"Cluster change summary saved to {summary_file}")
        else:
            logger.info("--summary requested but no history comparison was available — skipping")

    new_history_file = f'{args.output}_history.parquet.gz'
    new_db.to_parquet(new_history_file, compression='gzip')
    logger.info(f"History updated and saved to {new_history_file}")
    logger.info(f"Log saved to {log_file}")

if __name__=="__main__":
    main()
