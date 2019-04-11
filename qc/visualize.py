import io

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, average_precision_score


sns.set_context('paper')
sns.set_style('white')


def _output_figure(filename):
    out = None

    if filename is None:
        plt.show()
    elif filename == '__qcml_export__':
        binary = io.BytesIO()
        plt.savefig(binary, format='png')
        binary.seek(0)  # rewind the data
        out = binary.read()
    else:
        plt.savefig(filename)

    plt.close()

    return out


# Remember to use the Agg matplotlib backend for the heatmap annotation to work!
def plot_correlation_matrix(corr, filename=None):
    plt.figure(figsize=(11, 10))

    # generate a mask for the upper triangle
    mask = np.zeros_like(corr, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    sns.heatmap(corr, vmin=-1, vmax=1, linewidths=.5, square=True,
                xticklabels=corr.columns.values[:-1], yticklabels=corr.columns.values[1:],
                mask=mask, cbar_kws={'shrink': .75}, annot=True, fmt='.2f', annot_kws={'size': 4})

    # rotate overly long tick labels
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)

    return _output_figure(filename)


def _classes_to_colors(df):
    cmap = plt.cm.get_cmap('autumn')(np.linspace(0, 1, len(df.index.levels[0])))

    class_colors = {}
    color_idx = 0
    for c, _ in df.index.values:
        if class_colors.get(c) is None:
            class_colors[c] = cmap[color_idx]
            color_idx += 1

    colors = []
    for c, _ in df.index.values:
        colors.append(class_colors[c])

    return colors


def plot_timestamps(df, filename=None):
    plt.figure()

    plt.scatter(df.index.get_level_values(1), [0] * len(df.index.get_level_values(1)), 500, _classes_to_colors(df), '|')

    sns.despine(left=True)

    plt.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')

    return _output_figure(filename)


def _add_date_color_bar(df):
    num_ticks = 5
    ticker = mpl.ticker.MaxNLocator(num_ticks + 2, prune='both')

    mappable = cm.ScalarMappable(cmap=plt.cm.get_cmap('autumn'))
    mappable.set_array(range(num_ticks + 2))

    cb = plt.colorbar(mappable, ticks=ticker, shrink=0.75)
    cb.ax.set_yticklabels([df.index.values[i][1].strftime('%b %Y')
                           for i in range(0, len(df.index.values), len(df.index.values) // (num_ticks + 2))])
    cb.outline.set_linewidth(0)


def _scatter_plot(scatter_data, df, filename=None):
    plt.figure()

    plt.scatter(scatter_data[:, 0], scatter_data[:, 1], c=_classes_to_colors(df))

    sns.despine(left=True, bottom=True)

    plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    plt.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')

    _add_date_color_bar(df)

    return _output_figure(filename)


def plot_pca(df, filename=None):
    # transform data to lower dimension
    pca_data = PCA(2).fit_transform(df.values)

    # plot
    return _scatter_plot(pca_data, df, filename)


def plot_tsne(df, filename=None):
    # transform data to lower dimension
    tsne_data = TSNE(2, init='pca').fit_transform(df.values)

    # plot
    return _scatter_plot(tsne_data, df, filename)


def scatter_outliers(scatter_data, df, outlier_scores, score_threshold, filename=None):
    plt.figure()

    colors = _classes_to_colors(df)

    for i, d in enumerate(scatter_data):
        if outlier_scores[i] > score_threshold:
            size = 10 + (outlier_scores[i] - score_threshold) * 200
            color = colors[i]
            marker = 'o'
            alpha = None
        else:
            size = 10
            color = 'grey'
            marker = '.'
            alpha = 0.25
        plt.scatter(d[0], d[1], s=size, c=color, marker=marker, alpha=alpha)

    sns.despine(left=True, bottom=True)

    plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    plt.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')

    _add_date_color_bar(df)

    return _output_figure(filename)


def plot_pca_outliers(df, outlier_scores, score_threshold, filename=None):
    # transform data to lower dimension
    pca_data = PCA(2).fit_transform(df.values)

    # plot
    return scatter_outliers(pca_data, df, outlier_scores, score_threshold, filename)


def plot_tsne_outliers(df, outlier_scores, score_threshold, filename=None):
    # transform data to lower dimension
    tsne_data = TSNE(2, init='pca').fit_transform(df.values)

    # plot
    return scatter_outliers(tsne_data, df, outlier_scores, score_threshold, filename)


def plot_outlier_score_hist(outlier_scores, num_bins, score_cutoff, filename=None):
    plt.figure()

    ax = sns.distplot(outlier_scores, bins=np.arange(0, 1.01, 1 / num_bins), kde=False, axlabel='Outlier score (%)',
                      hist_kws={'histtype': 'stepfilled'})
    plt.ylabel('Number of experiments')

    if score_cutoff is not None:
        plt.axvline(score_cutoff, color=sns.color_palette()[0], linestyle='--')

    # convert axis labels to percentages
    ax.set_xticklabels([int(100 * float(label)) for label in ax.get_xticks()])

    sns.despine()

    return _output_figure(filename)


def plot_feature_importances(feature_importances, filename=None):
    sorted_importances = feature_importances.sort_values(ascending=False)

    with sns.axes_style('whitegrid'):
        fig = plt.figure()
        fig.set_tight_layout(True)

        sns.barplot(x=sorted_importances.index.values, y=sorted_importances, palette='Blues_d')

        plt.xticks(rotation='vertical', fontsize=5)

        return _output_figure(filename)


def plot_subspace_boxplots(data, highlights=None, filename=None):
    with sns.axes_style('whitegrid'):
        fig = plt.figure()
        fig.set_tight_layout(True)

        sns.boxplot(data=data, orient='v', palette='Blues_d')

        if highlights is not None:
            for i in range(len(highlights)):
                plt.plot(i, highlights[i], color='red', marker='d')

        plt.xticks(rotation=30, fontsize=10)

        return _output_figure(filename)


def plot_psm_boxplots(data, color_classes=None, filename=None, **kwargs):
    with sns.axes_style('whitegrid'):
        fig = plt.figure()
        fig.set_tight_layout(True)

        # add specific colors to the various box plots
        if color_classes is not None:
            kwargs['palette'] = [sns.color_palette()[c] for c in color_classes]

        sns.boxplot(data=data, **kwargs)

        if kwargs.get('orient') == 'h':
            plt.xlabel('Number of PSM\'s')
        else:
            plt.ylabel('Number of PSM\'s')

        return _output_figure(filename)


def plot_aucs(aucs, k_range, filename=None):
    max_auc = max(aucs)
    max_k = [k for k, a in zip(k_range, aucs) if a == max_auc]

    plt.figure()

    # plot all auc's
    plt.plot(k_range, aucs)
    # highlight max auc
    for k in max_k:
        plt.scatter(k, max_auc, s=50, c=sns.color_palette()[0], marker='D')

    plt.xlim(xmin=0)
    plt.ylim([0.5, 1.0])

    plt.xlabel('Local neighborhood size')
    plt.ylabel('AUC')

    return _output_figure(filename)


def plot_outlier_classes_score_hist(classes_scores, num_bins, filename=None):
    with sns.color_palette(sns.xkcd_palette(['medium green', 'orange yellow', 'faded red'])):
        plt.figure()

        # generate the histogram values for the three classes
        bins = np.arange(0, 1.01, 1 / num_bins)
        hist = pd.DataFrame({quality: np.histogram(classes_scores.loc[classes_scores['quality'] == quality]['score'],
                                                   bins=bins)[0] for quality in ['good', 'ok', 'poor']}, bins[:-1])

        ax = hist.plot(kind='bar', position=0)

        plt.xlabel('Outlier score (%)')
        plt.ylabel('Number of experiments')

        sns.despine(right=True, top=True)

        # change the x-axis to not include each bin value and convert to percentages
        ax.set_xticks(range(0, 21, 4))
        ax.set_xticklabels(range(0, 101, 20), rotation=0)

        return _output_figure(filename)


def plot_outlier_classes_score_kde(classes_scores, num_bins, filename=None):
    with sns.color_palette(sns.xkcd_palette(['medium green', 'orange yellow', 'faded red'])):
        plt.figure()

        bins = np.arange(0, 1.01, 1 / num_bins)
        for quality in ['good', 'ok', 'poor']:
            sns.distplot(classes_scores.loc[classes_scores['quality'] == quality]['score'], bins=bins, hist=False,
                         kde=True, kde_kws={'label': quality, 'shade': True}, norm_hist=True)

        plt.xlabel('Outlier score (%)')
        plt.ylabel('Density')

        sns.despine(right=True, top=True)

        # convert the outlier score tick labels to percentages
        ax = plt.gca()
        ax.set_xlim(0, 1)
        ax.set_xticklabels(range(0, 101, 20), rotation=0)

        return _output_figure(filename)


def _to_binary_class_labels(quality_classes, pos_label=('good', 'ok')):
    # convert quality classes: 1 -> poor, 0 -> good/ok
    # requires that NO unvalidated samples are present
    return np.array([0 if quality in pos_label else 1 for quality in quality_classes['quality']])


def plot_roc(classes_scores, filename=None):
    plt.figure()

    # convert ordinal class labels to binary labels
    binary_classes = _to_binary_class_labels(classes_scores)

    for zorder, ignore_quality in reversed(list(enumerate(['good', 'ok', None]))):
        # ignore samples of the quality that's not considered
        sample_weight = None if ignore_quality is None else _to_binary_class_labels(classes_scores, (ignore_quality, ))

        # compute roc
        fpr, tpr, _ = roc_curve(binary_classes, classes_scores['score'], sample_weight=sample_weight)
        auc = roc_auc_score(binary_classes, classes_scores['score'], sample_weight=sample_weight)
        # plot the ROC curve
        alpha = 1 if ignore_quality is None else 1 / 3
        label = 'all' if ignore_quality is None else "only '{}'".format('good' if ignore_quality == 'ok' else 'ok')
        plt.plot(fpr, tpr, zorder=zorder, alpha=alpha, label='ROC curve {} (AUC = {:.2f})'.format(label, auc))

    # plot the random ROC curve at 0.5
    plt.plot([0, 1], [0, 1], 'k--')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    plt.legend(loc='lower right')

    return _output_figure(filename)


def plot_precision_recall(classes_scores, filename=None):
    # compute false positive rate and true positive rate
    binary_classes = _to_binary_class_labels(classes_scores)
    precision, recall, _ = precision_recall_curve(binary_classes, classes_scores['score'])

    plt.figure()

    # plot the ROC curve
    plt.plot(recall, precision, label='Precision-recall curve (average precision = {:.2f})'
             .format(average_precision_score(binary_classes, classes_scores['score'])))

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])

    plt.xlabel('Recall')
    plt.ylabel('Precision')

    plt.legend(loc='lower right')

    return _output_figure(filename)


def plot_score_sensitivity_specificity(classes_scores, filename=None):
    # compute sensitivity and specificity
    sorted_scores = classes_scores.sort_values('score', ascending=False)
    sorted_binary_classes = _to_binary_class_labels(sorted_scores)

    sensitivity, specificity = np.zeros(len(sorted_scores)), np.zeros(len(sorted_scores))
    for i in range(len(sorted_binary_classes)):
        # true positives or false positives based on predictions above score cut-off
        tp = np.count_nonzero(sorted_binary_classes[:i + 1])
        fp = (i + 1) - tp
        # true negatives or false negatives based on predictions below score cut-off
        fn = np.count_nonzero(sorted_binary_classes[i + 1:])
        tn = len(sorted_binary_classes) - (i + 1) - fn

        # sensitivity and specificity
        sensitivity[i] = tp / (tp + fn)
        specificity[i] = tn / (tn + fp)

    plt.figure()
    ax1 = plt.gca()
    ax2 = plt.twinx()

    # plot the sensitivity and specificity in function of the outlier score
    p1 = ax1.plot(sorted_scores['score'], sensitivity, label='Sensitivity', color=sns.color_palette()[0])
    # advance colors for the second axis
    p2 = ax2.plot(sorted_scores['score'], specificity, label='Specificity', color=sns.color_palette()[1])

    ax1.set_xlim([-0.05, 1.05])
    ax1.set_ylim([-0.05, 1.05])
    ax2.set_ylim([-0.05, 1.05])

    ax1.set_xlabel('Outlier score')
    ax1.set_ylabel('Sensitivity')
    ax2.set_ylabel('Specificity')

    plots = p1 + p2
    ax1.legend(plots, [p.get_label() for p in plots], loc='center right')

    return _output_figure(filename)
