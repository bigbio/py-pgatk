import re
import pandas as pd
import datetime

from pypgatk.toolbox.general import ParameterConfiguration

class MzTabClassFdr(ParameterConfiguration):
    CONFIG_KEY_MzTabClassFdr = 'mzTab_class_fdr'
    CONFIG_DECOY_PREFIX = 'decoy_prefix'
    CONFIG_GLOBAL_FDR_CUTOFF = 'global_fdr_cutoff'
    CONFIG_CLASS_FDR_CUTOFF = 'class_fdr_cutoff'
    CONFIG_PEPTIDE_GROUPS_PREFIX = 'peptide_groups_prefix'

    def __init__(self, config_data, pipeline_arguments):
        """
      Init the class with the specific parameters.
      :param config_data configuration file
      :param pipeline_arguments pipelines arguments
      """

        super(MzTabClassFdr, self).__init__(self.CONFIG_KEY_MzTabClassFdr, config_data, pipeline_arguments)

        self._decoy_prefix = self.get_fdr_parameters(variable=self.CONFIG_DECOY_PREFIX, default_value='decoy')
        self._global_fdr_cutoff = self.get_fdr_parameters(variable=self.CONFIG_GLOBAL_FDR_CUTOFF, default_value=0.01)
        self._class_fdr_cutoff = self.get_fdr_parameters(variable=self.CONFIG_CLASS_FDR_CUTOFF, default_value=0.01)
        self._peptide_groups_prefix = self.get_fdr_parameters(variable=self.CONFIG_PEPTIDE_GROUPS_PREFIX,
                                                                        default_value={
                                                                            'non_canonical': ['altorf', 'pseudo',
                                                                                              'ncRNA'],
                                                                            'mutations': ['COSMIC', 'cbiomut'],
                                                                            'variants': ['var_mut', 'var_rs']})
        self._psm_search_engine_score_order =  {'1003113':True, '1003115': False,'1001493':True,'1001491':False}
                    
    def get_fdr_parameters(self, variable: str, default_value):
        value_return = default_value
        if variable in self.get_pipeline_parameters():
            value_return = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_MzTabClassFdr in self.get_default_parameters() and \
                variable in self.get_default_parameters()[self.CONFIG_KEY_MzTabClassFdr]:
            value_return = self.get_default_parameters()[self.CONFIG_KEY_MzTabClassFdr][variable]
        return value_return

    def _get_mzml_name(self, run, mtd):
        key = run + '-location'
        value = mtd.get(key)
        return value.split("/")[-1]

    def _is_decoy(self, accessions):
        list = accessions.split(',')
        if all(self._decoy_prefix in accession for accession in list):
            return 0
        else:
            return 1
    
    def _is_group(self, peptide_group_members, accessions):
        accession_group = 0
        list = accessions.split(',')
        for accession in list:
            for class_peptide in peptide_group_members:
                if class_peptide in accession:
                    accession_group += 1
        return len(list) == accession_group
    
    @staticmethod
    def _compute_global_fdr(df_psms, order):
        df_psms.sort_values("search_engine_score[1]", ascending=order, inplace=True)
        df_psms['FDR'] = (range(1, len(df_psms) + 1) / df_psms['target'].cumsum()) - 1
        df_psms['q-value'] = df_psms['FDR'][::-1].cummin()[::-1]

        df_psms.sort_values("search_engine_score[1]", ascending=order, inplace=True)

        return df_psms

    def _compute_class_fdr(self, df_psms,order):
        ls = []
        for c in self._peptide_groups_prefix:
            # split the dataframe and save the subset
            currClass = df_psms[
                df_psms['accession'].apply(lambda x: self._is_group(self._peptide_groups_prefix[c], x))]
            ls.append(currClass)

            # If there is no decoy to throw an exception
            if len(currClass[currClass["target"] == 0 ]) == 0:
                # raise ValueError(
                #     "There is not enough decoys to calculate " + c +" class-fdr.")
                print("Warning:There is no peptide or decoy of "+c+", and this kind of class-fdr has been skipped.")

            # calculate class-specific q-value
            currClass.sort_values("search_engine_score[1]", ascending=order, inplace=True)
            FDR = (range(1, len(currClass["target"]) + 1) / currClass["target"].cumsum()) - 1
            currClass['class-specific-q-value'] = FDR[::-1].cummin()[::-1]
        df = pd.concat(ls)

        # df_psms['class-specific-q-value'] = df['class-specific-q-value']
        df_psms = df_psms.merge(df['class-specific-q-value'], left_index=True, right_index=True, how='left')
        df_psms.loc[df_psms['class-specific-q-value'].isnull(), 'class-specific-q-value'] = df_psms['q-value']
        df_psms.sort_values("search_engine_score[1]", ascending=order, inplace=True)

        return df_psms

    def form_mzTab_class_fdr(self, input_mztab ,outfile_name):
        start_time = datetime.datetime.now()
        print("Start time :", start_time)
        
        file = open(input_mztab, "r")
        list = file.readlines()

        #Extract psms information
        psm = []
        mtd_dict = dict()
        for i in list:
            i = i.strip("\n")
            row_list = i.split('\t')
            if row_list[0] == "MTD":
                mtd_dict[row_list[1]] = row_list[2]
            elif row_list[0] == "PSH":
                psm_cols = row_list
            elif row_list[0] == "PSM":
                psm.append(row_list)

        psm_search_engine = mtd_dict.get("psm_search_engine_score[1]").split("MS:")[1][:7]
        order = self._psm_search_engine_score_order.get(psm_search_engine)

        #Convert to dataframe
        PSM = pd.DataFrame(psm, columns=psm_cols)
        PSM.loc[:, "SpecFile"] = PSM.apply(lambda x: self._get_mzml_name(x["spectra_ref"].split(":")[0], mtd_dict), axis=1)
        PSM.loc[:, "ScanNum"] = PSM.apply(lambda x: re.sub("[^\d]", "", x["spectra_ref"].split(":")[-1].split(" ")[-1]),axis=1)
        
        PSM.loc[:,"target"] = PSM.apply(lambda x: self._is_decoy(x["accession"]), axis=1)
        if len(PSM[PSM["target"] == 0]) ==0:
            raise ValueError(
                    "There is not enough decoys to calculate fdr.")

        PSM = self._compute_global_fdr(PSM, order)
        PSM = self._compute_class_fdr(PSM, order)
        PSM = PSM[((PSM['q-value'] < self._global_fdr_cutoff) & (
                PSM['class-specific-q-value'] < self._class_fdr_cutoff))]
        PSM.reset_index(drop=True, inplace=True)

        PSM.to_csv(outfile_name, header=1, sep="\t", index=None)
        
        end_time = datetime.datetime.now()
        print("End time :", end_time)
        time_taken = end_time - start_time
        print("Time consumption :", time_taken)




    