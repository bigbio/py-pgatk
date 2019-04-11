import math
import os
import urllib.request
import zipfile

import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.metrics import roc_auc_score

import export
import outlier
import qc_analysis
import visualize


##############################################################################

# OUTLIER VALIDATION BY PSM COMPARISON

def compare_outlier_psms(f_psms, outliers):
    # compare inliers and outliers based on their number of valid PSM's
    psms = pd.Series.from_csv(f_psms)

    outlier_psms = psms.filter(items=[index[0] for index in outliers.index.values])
    inlier_psms = psms.drop(outlier_psms.index)

    exporter.psm(inlier_psms, outlier_psms)

    return psms, inlier_psms, outlier_psms


def compare_outlier_subspace_psms(outliers, frequent_subspaces, psms, inlier_psms):
    # test whether a subspace can be related to a lower number of PSM's
    psm_table = pd.DataFrame(index=psms.index)
    psm_table['Inliers'] = inlier_psms
    color_classes = [0]
    pval_table = pd.DataFrame(index=range(len(frequent_subspaces)),
                              columns=['Metric(s)', 'Support (\%)', '\emph{p}-value'])
    for i, (subspace, support) in enumerate(frequent_subspaces):
        subspace = sorted(subspace)

        # compare outlier values
        outliers_values = pd.DataFrame([this_outlier for _, this_outlier in outliers.iterrows()
                                        if set(subspace) <= set(this_outlier.Subspace)])

        # compare outlier PSM's
        outlier_psms = psms.filter(items=[index[0] for index in outliers_values.index.values])

        # quantify difference between inliers and outliers
        t_stat, p_value = stats.ttest_ind(inlier_psms.values, outlier_psms.values, equal_var=False)

        psm_table['{}'.format(', '.join(subspace))] = outlier_psms
        color_classes.append(2 if p_value <= 0.05 and t_stat > 0 else 1)

        pval_table.set_value(i, 'Metric(s)', ', '.join(subspace))
        pval_table.set_value(i, 'Support (\%)', round(support))
        pval_table.set_value(i, '\emph{p}-value', '{}{:.5f}'.format('\cellcolor{lightgray} '
                                                                    if p_value <= 0.05 and t_stat > 0 else '', p_value))

    exporter.psm_pval(psm_table, pval_table, color_classes)


# OUTLIER VALIDATION USING PNNL EXPERT CLASSIFICATION


def find_optimal_outliers_k(data, f_class, k_min, dist):
    quality_classes = pd.DataFrame({'quality': pd.Series.from_csv(f_class)}, data.index.get_level_values(0))
    true_classes = visualize._to_binary_class_labels(quality_classes)
    k_range = np.arange(k_min, math.ceil(len(data) / 2), dtype=int)

    aucs = []
    for k in k_range:
        outlier_scores = outlier.detect_outliers_loop(data, k, metric=dist)
        aucs.append(roc_auc_score(true_classes, outlier_scores))
    max_auc = max(aucs)
    max_k = [k for k, auc in zip(k_range, aucs) if auc == max_auc]

    exporter.outlier_auc(aucs, k_range)

    return max_k, max_auc


def validate_outlier_score(data, f_class, scores, num_bins=20):
    # merge outlier scores and manual quality assignments
    classes_scores = pd.DataFrame({'quality': pd.Series.from_csv(f_class), 'score': scores})

    exporter.outlier_validation(classes_scores, num_bins)


##############################################################################


def run(args, f_psms=None, f_class=None, k_min=2, folder=None):
    global exporter
    exporter = export.Exporter(False, True, folder)
    qc_analysis.exporter = exporter

    data = qc_analysis.load_metrics(args.file_in, args.min_var, args.min_corr, args.scaling_mode)

    # compare outliers based on the number of psm's
    if f_psms is not None:
        outliers, outliers_score = qc_analysis.detect_outliers(data, args.k_neighbors, args.distance,
                                                               args.min_outlier, args.num_bins)
        frequent_subspaces = qc_analysis.analyze_outliers(data, outliers, args.k_neighbors, args.min_sup)
        psms, inlier_psms, outlier_psms = compare_outlier_psms(f_psms, outliers)
        compare_outlier_subspace_psms(outliers, frequent_subspaces, psms, inlier_psms)

    # compare outliers based on manual expert evaluation
    if f_class is not None:
        optimal_ks, _ = find_optimal_outliers_k(data, f_class, k_min, args.distance)
        outliers, outliers_score = qc_analysis.detect_outliers(data, optimal_ks[0], args.distance,
                                                               args.min_outlier, args.num_bins)
        validate_outlier_score(data, f_class, outliers_score, args.num_bins)
        qc_analysis.analyze_outliers(data, outliers, optimal_ks[0], args.min_sup)

    exporter.export(args.file_out)


if __name__ == '__main__':
    # check if the data folder exists and if all data files are present
    data_dir = 'data/'
    data_files = ['PNNL_iontrap_QuaMeter.tsv', 'PNNL_iontrap_validation.csv',
                  'PNNL_orbi_QuaMeter.tsv', 'PNNL_orbi_validation.csv',
                  'PNNL_velos_QuaMeter.tsv', 'PNNL_velos_validation.csv',
                  'TCGA_QuaMeter.tsv', 'TCGA_psms.csv']
    if not (os.path.exists(data_dir) and all([os.path.isfile(os.path.join(data_dir, f)) for f in data_files])):
        # otherwise download the data
        zip_filename = 'data.zip'
        urllib.request.urlretrieve('https://bitbucket.org/proteinspector/qc_analysis/downloads/data.zip', zip_filename)
        with zipfile.ZipFile(zip_filename) as zf:
            zf.extractall(data_dir)
        os.remove(zip_filename)

    # PNNL
    instruments = [('iontrap', 0.15), ('orbi', 0.25), ('velos', 0.20)]
    for instrument, outlier_score in instruments:
        run(qc_analysis.parse_args('-k 1 -o {} data/PNNL_{}_QuaMeter.tsv out.qcml'.format(outlier_score, instrument)),
            f_class='data/PNNL_{}_validation.csv'.format(instrument), folder='out/{}'.format(instrument))

    # TCGA
    run(qc_analysis.parse_args('-k 50 -o 0.25 data/TCGA_QuaMeter.tsv out.qcml'), f_psms='data/TCGA_psms.csv', folder='out/tcga')


##############################################################################
