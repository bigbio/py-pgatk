import argparse
import multiprocessing
import shlex

import fim
import pandas as pd




##############################################################################

# DATA LOADING AND PRE-PROCESSING
from qc import outlier, export


def load_metrics(file_in, min_var, min_corr, scaling):
    # load data from the input file
    from qc import preprocess
    data_raw = preprocess.load_metrics(file_in)

    # pre-process: remove low-variance and correlated metrics & scale the values
    data, variance, cor = preprocess.preprocess(data_raw, min_variance=min_var, min_corr=min_corr, scaling_mode=scaling)

    # add the preprocessing results to the qcML export
    exporter.low_variance(pd.Series(variance, index=data_raw.columns.values), min_var)
    exporter.correlation(cor, min_corr)
    exporter.preprocess_overview(data_raw.columns.values, variance, min_var, cor, min_corr)

    # add general visualizations to the qcML export
    exporter.global_visualization(data)

    return data


##############################################################################

# OUTLIER DETECTION

def detect_outliers(data, k, dist, outlier_threshold=None, num_bins=20):
    # compute outlier scores
    outlier_scores = outlier.detect_outliers_loop(data, k, metric=dist)

    # compute the outlier threshold (if required)
    if outlier_threshold is None:
        outlier_threshold = outlier.detect_outlier_score_threshold(outlier_scores, num_bins)

    # add the outlier score information to the qcML export
    exporter.outlier_scores(data, outlier_scores, outlier_threshold, num_bins)

    # separate significant outliers
    data_excluding_outliers, outliers = outlier.split_outliers(data, outlier_scores, outlier_threshold)

    return outliers, pd.Series(outlier_scores, data.index.get_level_values(0))


def analyze_outliers(data, outliers, k, min_sup):
    # retrieve explanatory subspaces for each outlier
    outliers['FeatureImportance'] = object
    outliers['Subspace'] = object
    with multiprocessing.Pool() as pool:
        # compute the subspace for each outlier
        subspaces = {name: pool.apply_async(outlier.get_outlier_subspace, (data, this_outlier, k))
                     for name, this_outlier in outliers.iterrows()}

        # set the outlier subspaces
        for name, result in subspaces.items():
            feature_importance, subspace = result.get()
            outliers.set_value(name, 'FeatureImportance', feature_importance.values)
            outliers.set_value(name, 'Subspace', subspace)

    # add the outliers' subspaces to the export
    for name, this_outlier in outliers.iterrows():
        exporter.outlier(this_outlier, data)

    # detect frequently occurring explanatory subspaces
    frequent_subspaces = sorted(fim.fim(outliers.Subspace, supp=min_sup, report='S'), key=lambda x: x[1], reverse=True)
    frequent_subspaces_table = pd.DataFrame(index=range(len(frequent_subspaces)),
                                            columns=['Outlier subspace QC metric(s)', 'Support (%)'])
    for i, (subspace, support) in enumerate(frequent_subspaces):
        frequent_subspaces_table.set_value(i, 'Outlier subspace QC metric(s)', ', '.join(subspace))
        frequent_subspaces_table.set_value(i, 'Support (%)', round(support))

    exporter.frequent_outlier_subspaces(frequent_subspaces_table, min_sup)

    return frequent_subspaces


##############################################################################

#  EXECUTE

# command-line execution
def parse_args(args_str=None):
    parser = argparse.ArgumentParser(description='Mass spectrometry quality control metrics analysis')
    parser.add_argument('file_in', type=argparse.FileType('r'),
                        help='the tab-separated input file containing the QC metrics')
    parser.add_argument('file_out', type=argparse.FileType('w'),
                        help='the name of the output file (.html extension for HTML export (default), '
                             '.qcml extension for qcML export')
    parser.add_argument('--min_var', '-var', default=0.0001, type=float,
                        help='metrics with a lower variance will be removed (default: %(default)s)')
    parser.add_argument('--min_corr', '-corr', default=0.9, type=float,
                        help='metrics with a higher correlation will be removed (default: %(default)s)')
    parser.add_argument('--scaling_mode', '-scale', default='robust', type=str, choices=['robust', 'standard'],
                        help='mode to standardize the metric values (default: %(default)s)')
    parser.add_argument('--k_neighbors', '-k', type=int, required=True,
                        help='the number of nearest neighbors used for outlier detection and interpretation')
    parser.add_argument('--distance', '-dist', default='manhattan', type=str,
                        help='metric to use for distance computation (default: %(default)s) '
                             'ny metric from scikit-learn or scipy.spatial.distance can be used')
    parser.add_argument('--min_outlier', '-o', default=None, type=float,
                        help='the minimum outlier score threshold (default: %(default)s) '
                             'if no threshold is provided, an automatic threshold is determined')
    parser.add_argument('--num_bins', '-bin', default=20, type=int,
                        help='the number of bins for the outlier score histogram (default: %(default)s)')
    parser.add_argument('--min_sup', '-sup', default=5, type=int,
                        help='the minimum support for subspace frequent itemset mining (default: %(default)s) '
                             'positive numbers are interpreted as percentages, negative numbers as absolute supports')

    # parse command-line arguments
    if args_str is None:
        return parser.parse_args()
    # or parse string arguments
    else:
        return parser.parse_args(shlex.split(args_str))


def run(args):
    global exporter
    exporter = export.Exporter(True, False)

    data = load_metrics(args.file_in, args.min_var, args.min_corr, args.scaling_mode)
    outliers, outliers_score = detect_outliers(data, args.k_neighbors, args.distance, args.min_outlier, args.num_bins)
    analyze_outliers(data, outliers, args.k_neighbors, args.min_sup)

    try:
        exporter.export(args.file_out)
    finally:
        args.file_out.close()


if __name__ == '__main__':
    run(parse_args())


##############################################################################
