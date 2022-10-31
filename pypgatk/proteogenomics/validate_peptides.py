import datetime
import os.path
import re
import pandas as pd
# from matplotlib import pyplot as plt
from pyopenms import *
from Bio import SeqIO

from pypgatk.toolbox.general import ParameterConfiguration

class ValidatePeptidesService(ParameterConfiguration):
    CONFIG_KEY_VALIDATE_PEPTIDES = 'validate_peptides'
    CONFIG_MZML_PATH = 'mzml_path'
    CONFIG_MZML_FILES = 'mzml_files'
    CONFIG_INFILE_NAME = 'infile_name'
    CONFIG_OUTFILE_NAME = 'outfile_name'
    CONFIG_IONS_TOLERANCE = 'ions_tolerance'
    CONFIG_RELATIVE = 'relative'
    CONFIG_MSGF = 'msgf'

    def __init__(self, config_data, pipeline_arguments):
        """
      Init the class with the specific parameters.
      :param config_data configuration file
      :param pipeline_arguments pipelines arguments
      """

        super(ValidatePeptidesService, self).__init__(self.CONFIG_KEY_VALIDATE_PEPTIDES, config_data, pipeline_arguments)

        self._ions_tolerance = self.get_validate_parameters(variable=self.CONFIG_IONS_TOLERANCE, default_value=0.02)
        self._relative = self.get_validate_parameters(variable=self.CONFIG_RELATIVE, default_value=False)
        self._msgf = self.get_validate_parameters(variable=self.CONFIG_MSGF, default_value=False)

    def get_validate_parameters(self, variable: str, default_value):
        value_return = default_value
        if variable in self.get_pipeline_parameters():
            value_return = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_VALIDATE_PEPTIDES in self.get_default_parameters() and \
                variable in self.get_default_parameters()[self.CONFIG_KEY_VALIDATE_PEPTIDES]:
            value_return = self.get_default_parameters()[self.CONFIG_KEY_VALIDATE_PEPTIDES][variable]
        return value_return

    # def _get_pep_pos(self, protein, sequence, fasta):
    #     for seq in SeqIO.parse(fasta, "fasta"):
    #         if seq.id == protein:
    #             sequence_canonical = str(seq.seq)
    #             rr = re.compile(sequence, re.I)
    #             match = re.finditer(rr, sequence_canonical)
    #             seq_position = int(re.sub("\D", "", protein.split(":")[2]))
    #             for i in match:
    #                 start = i.start() + 1
    #                 end = i.end()
    #                 if seq_position >= start and seq_position <= end:
    #                     return seq_position - start + 1
    #                 else:
    #                     continue

    # def _get_canonical_aa(self, protein):
    #     return protein.split(":")[2].split(".")[1][0]

    # def _get_variant_aa(self, protein):
    #     return protein.split(":")[2].split(".")[1][-1]

    # def _get_canonical_peptide(self, variant_peptide,canonical_aa,position):
    #     if position==0:
    #         return "Unmutated"
    #     else:
    #         str1 = variant_peptide[0:position - 1]
    #         str2 = variant_peptide[position:]
    #         return str1 + canonical_aa + str2

    # def get_position(self, input_psm_table, input_fasta, output_psm_table):
    #     if self._msgf:
    #         PSM = pd.read_table(input_psm_table, header = 0, sep = "\t", index_col=0)
    #         PSM.loc[:, "Variant Peptide"] = PSM.apply(lambda x: re.sub("[^A-Z]","",x["Peptide"]), axis = 1)
    #         PSM.loc[:, "Canonical AA"] = PSM.apply(lambda x: self._get_canonical_aa(x["Protein"]), axis = 1)
    #         PSM.loc[:, "Variant AA"] = PSM.apply(lambda x: self._get_variant_aa(x["Protein"]), axis = 1)
    #         PSM.loc[:, "position"] = PSM.apply(lambda x: self._get_pep_pos(x["Protein"], x["Variant Peptide"], input_fasta), axis = 1)
    #         PSM["position"].fillna(0, inplace = True)

    #         PSM.loc[:, "Canonical Peptide"] = PSM.apply(lambda x: self._get_canonical_peptide(x["Variant Peptide"],x["Canonical AA"],int(x["position"])), axis = 1)
    #         psm = list(PSM)
    #         PSM = PSM.loc[:,psm]
    #     else:
    #         PSM = pd.read_table(input_psm_table, header = 0, sep = "\t", index_col=0)
    #         PSM.loc[:, "Canonical AA"] = PSM.apply(lambda x: self._get_canonical_aa(x["accession"]), axis = 1)
    #         PSM.loc[:, "Variant AA"] = PSM.apply(lambda x: self._get_variant_aa(x["accession"]), axis = 1)
    #         PSM.loc[:, "position"] = PSM.apply(lambda x: self._get_pep_pos(x["accession"], x["sequence"], input_fasta), axis = 1)
    #         PSM["position"].fillna(0, inplace = True)

    #     PSM.to_csv(output_psm_table, sep = "\t", index = 0)

    def _predict_MS2_spectrum(self, Peptide, size, product_ion_charge = 1 ):
        if self._msgf:
            Peptide = re.sub("[-?]", "", Peptide)
            modification = re.finditer("(\+\d{1,}\.\d{1,})", Peptide)

            a = 0
            for i in modification:
                Peptide = Peptide[:i.start() + a] + '[' + Peptide[i.start() + a:i.end() + a] + ']' + Peptide[i.end() + a:]
                a += 2

        tsg = TheoreticalSpectrumGenerator()
        spec = MSSpectrum()
        peptide = AASequence.fromString(Peptide)
        # size = len(peptide.toUnmodifiedString())
        p = Param()
        p.setValue("add_metainfo", "true")
        p.setValue("add_first_prefix_ion", "true")
        p.setValue("add_precursor_peaks", "true")
        tsg.setParameters(p)
        tsg.getSpectrum(spec, peptide, 1, 1)  # charge range 1:1

        b_y_ions = []
        for i in spec.getStringDataArrays()[0]:
            b_y_ions.append(i.decode())
        mz = []
        for i in spec:
            mz.append(i.getMZ())

        ions = pd.DataFrame({"mz":mz,"ion":b_y_ions,"z":1})

        ions.loc[2 * size - 2, "ion"] = "b" + str(size)
        ions = ions.drop(2 * size - 1)
        ions.loc[2 * size, "ion"] = "y" + str(size)

        ions.loc[:, "ion"] = ions.apply(lambda x: re.sub("[+]", "", x["ion"]), axis=1)
        ions.loc[:, "pos"] = ions.apply(lambda x: re.sub("[^\d]", "", x["ion"]), axis=1)
        ions.loc[:, "type"] = ions.apply(lambda x: re.sub("[^a-z]", "", x["ion"]), axis=1)

        proton_mono_mass = 1.007276
        if product_ion_charge > 1:
            ions2 = ions.copy()
            ions2.loc[:, "mz"] = ions2.apply(lambda x: (x["mz"]+proton_mono_mass)/2, axis=1)
            ions2.loc[:, "z"] = 2

            ions = ions.merge(ions2, how='outer')

        ions = ions.reset_index(drop=True)

        return ions    

    def _get_intensity(self, exp_peak, ion_mz):
        exp_peak.loc[:,"mz_difference"] = exp_peak.apply(lambda x:abs(float(ion_mz) - x["mz"]), axis = 1)
        min_index=exp_peak["mz_difference"].idxmin()
        return exp_peak.loc[exp_peak["mz_difference"]==exp_peak["mz_difference"].min()].loc[min_index,"intensity"]

    def _match_exp2predicted(self, exp_peak, pred_peak):
        pred_peak.loc[:,"error"] = pred_peak.apply(lambda x:min(abs(float(x["mz"])-exp_peak["mz"])), axis = 1)
        pred_peak.loc[:,"intensity"] = pred_peak.apply(lambda x:self._get_intensity(exp_peak,x["mz"]), axis = 1)
        pred_peak.loc[:,"ppm"] = pred_peak.apply(lambda x:round(x["error"]/x["mz"]*1000000,2), axis = 1)

        if self._relative:
            match_ions = pred_peak[pred_peak["ppm"] < self._ions_tolerance]
        else:
            match_ions = pred_peak[pred_peak["error"] < self._ions_tolerance]

        match_ions = match_ions.reset_index(drop=True)

        return match_ions

    def _InspectSpectrum(self, DF, mzml_path, mzml_files):
        if self._msgf:
            DF.loc[:, "peptide_length"] = DF.apply(lambda x: len(re.sub("[^A-Z]", "", x["Peptide"])), axis=1)
        else:
            DF.loc[:, "peptide_length"] = DF.apply(lambda x: len(x["sequence"]), axis=1)

        DF["status"] = "skiped"

        DF["ions_support"] = "NO"
        DF["support_ions"] = ""
        DF["sum.supportions.intensity"] = float(0)

        DF["flanking_ions_support"] = "NO"
        DF["flanking_ions"] = ""
        DF["sum.flanking.ions.intensity"] = float(0)

        DF["matched_ions"] = ""
        DF["sum.matchedions.intensity"] = float(0)
        DF["sum.fragmentions.intensity"] = float(0)
        DF["maxintensity"] = float(0)
        DF["average_intensity"] = float(0)
        DF["median_intensity"] = float(0)

        Spectra_list = {}
        for k in range(DF["SpecFile"].nunique()):
            Spectra_list[DF["SpecFile"].unique()[k]] = []

        for i in range(DF.shape[0]):
            spectra_file = str(DF.loc[i, "SpecFile"])

            if mzml_files and not mzml_path:
                mzml_list = mzml_files.split(",")
                for file in mzml_list:
                    if spectra_file in file:
                        mzml_file = file
                        break
            elif not mzml_files and mzml_path:
                mzml_file = os.path.join(mzml_path, spectra_file)
            else:
                raise ValueError(
                    "You only need to use either '--mzml_path' or '--mzml_files'.")

            scan_num = int(DF.loc[i, "ScanNum"])
            if self._msgf:
                # seq = DF.loc[i, "Variant Peptide"]
                seq = re.sub("[^A-Z]", "", DF.loc[i, "Peptide"])
                length = DF.loc[i, "peptide_length"]
            else:
                seq = DF.loc[i, "sequence"]
                length = DF.loc[i, "peptide_length"]

            if not Spectra_list[spectra_file]:
                # read the mass spectrometry file
                exp = MSExperiment()
                MzMLFile().load(mzml_file, exp)
                Spectra_list[spectra_file].append(exp)
                look = SpectrumLookup()
                look.readSpectra(exp, "((?<SCAN>)\d+$)")
                Spectra_list[spectra_file].append(look)
            # get peaks through ScanNum
            index = Spectra_list[spectra_file][1].findByScanNumber(scan_num)
            exp_peaks = Spectra_list[spectra_file][0].getSpectrum(index).get_peaks()
            exp_peaks = pd.DataFrame({"mz": exp_peaks[0], "intensity": exp_peaks[1]})

            if self._msgf:
                predicted_peaks = self._predict_MS2_spectrum(str(DF.loc[i, "Peptide"]), length, 1)
            else:
                predicted_peaks = self._predict_MS2_spectrum(str(DF.loc[i, "opt_global_cv_MS:1000889_peptidoform_sequence"]), length, 1)
            match_ions = self._match_exp2predicted(exp_peaks, predicted_peaks)

            maxintensity = exp_peaks["intensity"].max()
            average_intensity = exp_peaks["intensity"].mean()
            median_intensity = exp_peaks["intensity"].median()

            DF.loc[i, "sum.fragmentions.intensity"] = exp_peaks["intensity"].sum()
            DF.loc[i, "maxintensity"] = maxintensity
            DF.loc[i, "average_intensity"] = average_intensity
            DF.loc[i, "median_intensity"] = median_intensity

            if match_ions.shape[0] == 0:
                continue
            DF.loc[i, "matched_ions"] = ','.join(match_ions["ion"].unique().tolist())
            DF.loc[i, "sum.matchedions.intensity"] = match_ions["intensity"].sum()

            if DF.loc[i, "position"] == "canonical":
                continue
            if DF.loc[i, "position"] == "non-canonical":
                continue
            position = int(DF.loc[i, "position"])
            if position == 0:
                continue
            if position > length:
                continue
            
            DF.loc[i, "status"] = "checked"
            supportions_intensity = 0
            ions_support = "NO"
            supportions = ""

            for j in range(match_ions.shape[0]):
                type = match_ions.loc[j, "type"]
                pos = int(match_ions.loc[j, "pos"])
                ion = match_ions.loc[j, "ion"]

                if type == "b" and pos >= position:
                    ions_support = "YES"
                    supportions_intensity = supportions_intensity + match_ions.loc[j, "intensity"]
                    supportions = supportions + ',' + ion
                elif type == "y" and pos > length - position:
                    ions_support = "YES"
                    supportions_intensity = supportions_intensity + match_ions.loc[j, "intensity"]
                    supportions = supportions + ',' + ion

            DF.loc[i, "ions_support"] = ions_support
            DF.loc[i, "support_ions"] = supportions
            DF.loc[i, "sum.supportions.intensity"] = supportions_intensity

            # check if it is a noise peak or isotope peak supporting mutant ions
            if DF.loc[i, "sum.supportions.intensity"] < DF.loc[i, "median_intensity"]:
                DF.loc[i, "ions_support"] = "NO"

            flanking_ions_support = "NO"
            n1 = length
            n2 = position
            match_ions_set = set(match_ions["ion"].tolist())

            if n2 == 1:
                flanking_ions = {"b1", "y" + str(n1 - 1)}
                flanking_ions = flanking_ions.intersection(match_ions_set)
                if len(flanking_ions) > 0:
                    flanking_ions_support = "YES"
            elif n2 == n1:
                flanking_ions = {"y1", "b" + str(n1 - 1)}
                flanking_ions = flanking_ions.intersection(match_ions_set)
                if len(flanking_ions) > 0:
                    flanking_ions_support = "YES"
            else:
                flanking_ions_left = {"b" + str(n2 - 1), "y" + str(n1 - n2 + 1)}
                flanking_ions_right = {"b" + str(n2), "y" + str(n1 - n2)}

                flanking_ions_left = flanking_ions_left.intersection(match_ions_set)
                flanking_ions_right = flanking_ions_right.intersection(match_ions_set)

                flanking_ions = flanking_ions_left.union(flanking_ions_right)
                if len(flanking_ions_left) > 0 and len(flanking_ions_right) > 0:
                    flanking_ions_support = "YES"

            DF.loc[i, "flanking_ions_support"] = flanking_ions_support
            DF.loc[i, "flanking_ions"] = ",".join(flanking_ions)
            if flanking_ions:
                DF.loc[i, "sum.flanking.ions.intensity"] = match_ions[match_ions['ion'].str.contains("|".join(flanking_ions))]["intensity"].sum()

            if DF.loc[i, "sum.flanking.ions.intensity"] < DF.loc[i, "median_intensity"]:
                DF.loc[i, "flanking_ions_support"] = "NO"

            # fragmentation is not preferable at Cterm side of proline, so only require supporting ions
            if re.search("P", seq[position - 1:position]):
                DF.loc[i, "flanking_ions_support"] = DF.loc[i, "ions_support"]

        return DF

    def validate(self, infile_name, outfile_name, mzml_path, mzml_files):
        start_time = datetime.datetime.now()
        print("Start time :", start_time)
        df_psm = pd.read_table(infile_name, header=0, sep="\t")
        df_output = self._InspectSpectrum(df_psm, mzml_path, mzml_files)
        df_output.to_csv(outfile_name, header=1, sep="\t")

        # if self._msgf:
        #     df_sub = df_output[df_output["status"] == "checked"]
        #     saav_psm_passed = df_sub[df_sub["flanking_ions_support"]=="YES"]["PrecursorError(ppm)"]
        #     saav_psm_failed = df_sub[df_sub["flanking_ions_support"]=="NO"]["PrecursorError(ppm)"]
        #     plot=plt.figure(figsize=(10,7))
        #     plot1=plot.add_subplot(1,2,1)
        #     plot2=plot.add_subplot(1,2,2)
        #     plot1.hist(saav_psm_passed,bins=20)
        #     plot1.set_xlabel("PrecursorError(ppm)")
        #     plot1.set_title("SpectrumAI curated")
        #     plot2.hist(saav_psm_failed,bins=20)
        #     plot2.set_xlabel("PrecursorError(ppm)")
        #     plot2.set_title("SpectrumAI discarded")
        #     plt.savefig("precursorError_histogram.pdf")

        end_time = datetime.datetime.now()
        print("End time :", end_time)
        time_taken = end_time - start_time
        print("Time consumption :", time_taken)



        