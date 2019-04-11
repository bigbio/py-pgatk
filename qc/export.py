import base64
import datetime
import itertools
import os
import sqlite3

import jinja2
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from qc import qcml, visualize


def extract_idp_psms_to_file(f_in, f_out):
    # export the number of PSM's for each sample from the IDPicker database to a csv file
    conn = sqlite3.connect(f_in)

    psms = {}
    c = conn.cursor()
    for result in c.execute(
            'SELECT SS.Name, COUNT(*) FROM PeptideSpectrumMatch PSM, Spectrum S, SpectrumSource SS WHERE PSM.Spectrum = S.Id AND S.Source = SS.Id GROUP BY SS.Id'):
        psms[result[0]] = result[1]

    pd.Series(psms).to_csv(f_out)


class Exporter:

    def __init__(self, export_report=True, export_figures=False, fig_folder=None):
        self.export_report = export_report
        self.export_figures = export_figures
        self.fig_folder = fig_folder

        self.creation_date = datetime.datetime.now()

        if export_figures and self.fig_folder is not None:
            if not os.path.exists(self.fig_folder):
                os.makedirs(self.fig_folder)
        elif export_figures:
            self.fig_folder = os.getcwd()

        if self.export_report:
            # nothing needed to be done for HTML export, we will differ between HTML and qcML later on

            # create qcML
            self.qcml_out = qcml.qcMLType()
            self.qcml_out.set_version('0.0.8')

            # add references to the controlled vocabularies (CV's)
            self.cv_ms = qcml.CVType('PSI-MS', '3.78.0', 'http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo', 'MS')
            self.cv_qc = qcml.CVType('MS-QC', '0.1.1', 'https://github.com/qcML/qcML-development/blob/master/cv/qc-cv.obo', 'QC')
            self.cv_outlier = qcml.CVType('Outlier detection and interpretation', '0', 'article/code link', 'outlier')

            self.qcml_out.set_cvList(qcml.CVListType([self.cv_ms, self.cv_qc, self.cv_outlier]))

            # add global outlier information in a setQuality
            self.set_quality = qcml.SetQualityAssessmentType(ID='OutlierAnalysis')
            self.set_quality.add_metaDataParameter(qcml.MetaDataType(name='Creation date', value=self.creation_date,
                                                                     cvRef=self.cv_qc.get_ID(), accession='MS:1000747',
                                                                     ID='{}_CreationDate'.format(
                                                                         self.set_quality.get_ID())))
            self.qcml_out.add_setQuality(self.set_quality)

            # add embedded stylesheet
            with open('to_html.xsl', 'r') as xsl:
                self.qcml_out.set_embeddedStylesheetList(qcml.embeddedStylesheetListType(anytypeobjs_=xsl.read()))

        if self.export_figures:
            pass

    def low_variance(self, variances, min_var):
        if self.export_report:
            param_var = qcml.QualityParameterType(name='Variance threshold', ID='VarianceThreshold', value='{:.3e}'.format(min_var),
                                                  cvRef=self.cv_outlier.get_ID(), accession='none')
            self.set_quality.add_qualityParameter(param_var)

            values = [' '.join((v[0], '{:.3e}'.format(v[1]))) for v in variances[variances <= min_var].iteritems()]
            table = qcml.TableType(tableColumnTypes=['RemovedMetric', 'Variance'], tableRowValues=values)
            self.set_quality.add_attachment(qcml.AttachmentType(name='Low variance metrics', ID='var', table=table,
                                                                cvRef=self.cv_outlier.get_ID(), accession='none',
                                                                qualityParameterRef=param_var.get_ID()))
        if self.export_figures:
            pass

    def correlation(self, corr, min_corr):
        if self.export_report:
            param_corr = qcml.QualityParameterType(name='Correlation threshold', ID='CorrelationThreshold', value=min_corr,
                                                   cvRef=self.cv_outlier.get_ID(), accession='none')
            self.set_quality.add_qualityParameter(param_corr)

            values = []
            corr_features = set()
            for row in range(len(corr.index)):
                if corr.columns.values[row] not in corr_features:
                    for col in range(row + 1, len(corr.columns)):
                        if corr.columns.values[col] not in corr_features and abs(corr.iloc[row, col]) > min_corr:
                            corr_features.add(corr.columns.values[col])
                            values.append(' '.join((corr.columns.values[row], corr.columns.values[col], '{:.2%}'.format(corr.iloc[row, col]))))
            table = qcml.TableType(tableColumnTypes=['RetainedMetric', 'RemovedMetric', 'Correlation'], tableRowValues=values)
            self.set_quality.add_attachment(qcml.AttachmentType(name='Correlated metrics', ID='corr', table=table,
                                                                cvRef=self.cv_outlier.get_ID(), accession='none',
                                                                qualityParameterRef=param_corr.get_ID()))

        if self.export_figures:
            visualize.plot_correlation_matrix(corr, os.path.join(self.fig_folder, 'corr.pdf'))

    def preprocess_overview(self, metrics, variances, min_var, correlation, min_corr):
        if self.export_figures:
            with open(os.path.join(self.fig_folder, 'table_preprocess.txt'), 'w') as f_out:
                table = pd.DataFrame(index=range(len(metrics)), columns=['Metric', '{Variance}', 'Correlated Metric', '{Correlation (\%)}'])

                table['Metric'] = metrics

                table['{Variance}'] = variances

                corr_features = set()
                for row in range(len(correlation.index)):
                    if correlation.columns.values[row] not in corr_features:
                        for col in range(row + 1, len(correlation.columns)):
                            if correlation.columns.values[col] not in corr_features and abs(correlation.iloc[row, col]) > min_corr:
                                retained_metric = correlation.columns.values[row]
                                removed_metric = correlation.columns.values[col]
                                removed_idx = np.where(metrics == removed_metric)[0][0]
                                corr_features.add(removed_metric)
                                table.ix[removed_idx, '{Correlation (\%)}'] = '{:.4}'.format(correlation.iloc[row, col] * 100)
                                table.ix[removed_idx, 'Correlated Metric'] = retained_metric

                f_out.write(table.to_latex(index=False, escape=False, na_rep='',
                                           formatters={'{Variance}': lambda x: '{}{:.2e}'.format('\cellcolor{gray} ' if x < min_var else '', x)}))

    def global_visualization(self, data):
        if self.export_report:
            self.set_quality.add_attachment(qcml.AttachmentType(name='Experiment execution time', ID='time',
                                                                binary=visualize.plot_timestamps(data, filename='__qcml_export__'),
                                                                cvRef=self.cv_outlier.get_ID(), accession='none'))
            self.set_quality.add_attachment(qcml.AttachmentType(name='PCA visualization', ID='PCA',
                                                                binary=visualize.plot_pca(data, filename='__qcml_export__'),
                                                                cvRef=self.cv_outlier.get_ID(), accession='none'))
            self.set_quality.add_attachment(qcml.AttachmentType(name='t-SNE visualization', ID='t-SNE',
                                                                binary=visualize.plot_tsne(data, filename='__qcml_export__'),
                                                                cvRef=self.cv_outlier.get_ID(), accession='none'))

        if self.export_figures:
            visualize.plot_timestamps(data, filename=os.path.join(self.fig_folder, 'dates.pdf'))
            visualize.plot_pca(data, filename=os.path.join(self.fig_folder, 'pca.pdf'))
            visualize.plot_tsne(data, filename=os.path.join(self.fig_folder, 'tsne.pdf'))

            pca = PCA(2)
            pca.fit_transform(data.values)

            with open(os.path.join(self.fig_folder, 'table_pca.txt'), 'w') as f_out:
                pca_table = pd.DataFrame(index=range(len(data.columns.values)), columns=['Metric', 0, 1])

                pca_table['Metric'] = data.columns.values
                pca_table[0] = pca.components_[0]
                pca_table[1] = pca.components_[1]
                pca_table.columns = ['Metric',
                                     '{{PC 1 ({:.1f}\,\%)}}'.format(pca.explained_variance_ratio_[0] * 100),
                                     '{{PC 2 ({:.1f}\,\%)}}'.format(pca.explained_variance_ratio_[1] * 100)]

                f_out.write(pca_table.to_latex(index=False, escape=False, float_format=lambda x: '{}{:.5f}'.format(
                    '\cellcolor{gray} ' if abs(x) >= 0.3 else '\cellcolor{lightgray} ' if abs(x) >= 0.2 else '', x)))

    def outlier_scores(self, data, outlier_scores, outlier_threshold, num_bins):
        if self.export_report:
            param_score = qcml.QualityParameterType(name='Outlier score threshold', ID='OutlierScoreThreshold', value=outlier_threshold,
                                                    cvRef=self.cv_outlier.get_ID(), accession='none')
            self.set_quality.add_qualityParameter(param_score)
            param_nr = qcml.QualityParameterType(name='Number of outliers', ID='NrOutliers', value=(outlier_scores > outlier_threshold).sum(),
                                                 cvRef=self.cv_outlier.get_ID(), accession='none')
            self.set_quality.add_qualityParameter(param_nr)

            attach_hist = qcml.AttachmentType(name='Outlier score histogram', ID='OutlierScoreHistogram', qualityParameterRef=param_score.get_ID(),
                                              binary=visualize.plot_outlier_score_hist(outlier_scores, num_bins, outlier_threshold, filename='__qcml_export__'),
                                              cvRef=self.cv_outlier.get_ID(), accession='none')
            self.set_quality.add_attachment(attach_hist)

        if self.export_figures:
            visualize.plot_outlier_score_hist(outlier_scores, num_bins, outlier_threshold, filename=os.path.join(self.fig_folder, 'outlier-hist.pdf'))
            visualize.plot_pca_outliers(data, outlier_scores, outlier_threshold, filename=os.path.join(self.fig_folder, 'pca-outlier.pdf'))
            visualize.plot_tsne_outliers(data, outlier_scores, outlier_threshold, filename=os.path.join(self.fig_folder, 'tsne-outlier.pdf'))

    def outlier(self, outlier, data):
        feature_importance = pd.Series(outlier['FeatureImportance'], index=outlier.drop(['OutlierScore', 'FeatureImportance', 'Subspace']).index)

        if self.export_report:
            run_quality = qcml.RunQualityAssessmentType(ID=outlier.name[0])
            self.qcml_out.add_runQuality(run_quality)

            run_quality.add_metaDataParameter(qcml.MetaDataType(name='Creation date', value=self.creation_date,
                                                                cvRef=self.cv_qc.get_ID(), accession='MS:1000747',
                                                                ID='{}_CreationDate'.format(run_quality.get_ID())))

            score = qcml.QualityParameterType(name='Outlier score', value=outlier['OutlierScore'],
                                              ID='{}_OutlierScore'.format(run_quality.get_ID()),
                                              cvRef=self.cv_outlier.get_ID(), accession='none')
            run_quality.add_qualityParameter(score)

            fig_features = visualize.plot_feature_importances(feature_importance, filename='__qcml_export__')
            fig_subspace = visualize.plot_subspace_boxplots(data[outlier['Subspace']], outlier[outlier['Subspace']], filename='__qcml_export__')

            run_quality.add_attachment(qcml.AttachmentType(name='Feature importance', binary=fig_features,
                                                           ID='{}_FeatureImportance'.format(run_quality.get_ID()),
                                                           cvRef=self.cv_outlier.get_ID(), accession='none',
                                                           qualityParameterRef=score.get_ID()))
            run_quality.add_attachment(qcml.AttachmentType(name='Explanatory subspace', binary=fig_subspace,
                                                           ID='{}_Subspace'.format(run_quality.get_ID()),
                                                           cvRef=self.cv_outlier.get_ID(), accession='none',
                                                           qualityParameterRef=score.get_ID()))

        if self.export_figures:
            if not os.path.exists(os.path.join(self.fig_folder, 'outlier/')):
                os.makedirs(os.path.join(self.fig_folder, 'outlier/'))
            visualize.plot_feature_importances(feature_importance, filename=os.path.join(self.fig_folder, 'outlier/{}_features.pdf'.format(outlier.name[0])))
            visualize.plot_subspace_boxplots(data[outlier['Subspace']], outlier[outlier['Subspace']], filename=os.path.join(self.fig_folder, 'outlier/{}_subspace.pdf'.format(outlier.name[0])))

    def frequent_outlier_subspaces(self, subspaces, min_sup):
        if self.export_report:
            if min_sup > 0:
                min_sup_str = '{}%'.format(min_sup)
            else:
                min_sup_str = '{}'.format(min_sup * -1)
            param_support = qcml.QualityParameterType(name='Minimum support', ID='minsup',
                                                      value=min_sup_str,
                                                      cvRef=self.cv_outlier.get_ID(), accession='none')
            self.set_quality.add_qualityParameter(param_support)

            values = ['{} {}'.format(subspace.iloc[0].replace(', ', '_'), subspace.iloc[1]) for _, subspace in subspaces.iterrows()]
            table = qcml.TableType(tableColumnTypes=['Subspace', 'Support(%)'], tableRowValues=values)
            self.set_quality.add_attachment(qcml.AttachmentType(name='Frequently occuring explanatory subspaces', ID='freq',
                                                                table=table, qualityParameterRef=param_support.get_ID(),
                                                                cvRef=self.cv_outlier.get_ID(), accession='none'))

        if self.export_figures:
            with open(os.path.join(self.fig_folder, 'table_freq.txt'), 'w') as f_out:
                f_out.write(subspaces.to_latex(index=False))

    def psm(self, inlier_psms, outlier_psms):
        if self.export_report:
            pass

        if self.export_figures:
            visualize.plot_psm_boxplots(pd.DataFrame({'Inliers ({})'.format(len(inlier_psms)): inlier_psms,
                                                      'Outliers ({})'.format(len(outlier_psms)): outlier_psms}),
                                        filename=os.path.join(self.fig_folder, 'psm_all.pdf'))

    def psm_pval(self, psms, pvals, color_classes):
        if self.export_report:
            pass

        if self.export_figures:
            with open(os.path.join(self.fig_folder, 'table_psm_pval.txt'), 'w') as f_out:
                f_out.write(pvals.to_latex(index=False, escape=False))

            visualize.plot_psm_boxplots(psms, color_classes, orient='h', filename=os.path.join(self.fig_folder, 'psm_subspace.pdf'))

    def outlier_auc(self, aucs, k_range):
        if self.export_figures:
            visualize.plot_aucs(aucs, k_range, filename=os.path.join(self.fig_folder, 'auc.pdf'))

    def outlier_validation(self, classes_scores, num_bins):
        if self.export_figures:
            visualize.plot_outlier_classes_score_hist(classes_scores, num_bins, os.path.join(self.fig_folder, 'outlier-scorequality-hist.pdf'))
            visualize.plot_outlier_classes_score_kde(classes_scores, num_bins, os.path.join(self.fig_folder, 'outlier-scorequality-kde.pdf'))
            visualize.plot_roc(classes_scores, os.path.join(self.fig_folder, 'roc.pdf'))
            visualize.plot_precision_recall(classes_scores, os.path.join(self.fig_folder, 'precision-recall.pdf'))
            visualize.plot_score_sensitivity_specificity(classes_scores, os.path.join(self.fig_folder, 'score-sensitivity-specificity.pdf'))

    def export(self, file_out):
        if self.export_report:
            _, ext = os.path.splitext(file_out.name)
            if ext == 'html':
                self.export_to_html(file_out)
            elif ext == 'qcml':
                self.qcml_out.export(file_out, 0, name_='qcML', namespacedef_='xmlns="http://www.prime-xs.eu/ms/qcml"')
            else:       # default: export to HTML
                self.export_to_html(file_out)

        if self.export_figures:
            pass

    def export_to_html(self, file_out):
        context = {}
        # summary information
        for parameter in itertools.chain(self.set_quality.get_qualityParameter(), self.set_quality.get_attachment()):
            identifier = ''.join(ch for ch in parameter.ID if ch.isalnum()).lower()
            param_dict = {'name': parameter.name, 'value': parameter.value}

            if hasattr(parameter, 'binary') and parameter.binary is not None:
                param_dict['binary'] = base64.b64encode(parameter.binary).decode('ascii')
            elif hasattr(parameter, 'table') and parameter.table is not None:
                param_dict['header'] = parameter.table.tableColumnTypes
                param_dict['rows'] = parameter.table.tableRowValues

            context[identifier] = param_dict

        # individual outliers
        outliers = []
        for outlier in self.qcml_out.get_runQuality():
            outlier_dict = {'name': outlier.ID, 'score': outlier.get_qualityParameter()[0].value}

            for attachment in outlier.get_attachment():
                if attachment.name == 'Feature importance':
                    outlier_dict['features'] = base64.b64encode(attachment.binary).decode('ascii')
                elif attachment.name == 'Explanatory subspace':
                    outlier_dict['subspace'] = base64.b64encode(attachment.binary).decode('ascii')

            outliers.append(outlier_dict)

        context['outliers'] = outliers

        file_out.write(jinja2.Environment(loader=jinja2.FileSystemLoader('.')).get_template('template.html').render(context))
