import os

import gffutils
import vcf
from Bio import SeqIO
from Bio.Seq import Seq

from db.VCFtoProteinDB import get_features, check_overlap, get_altseq, get_orfs
from toolbox.general import ParameterConfiguration


class EnsemblDataService(ParameterConfiguration):
    NUCLEAR_TRANSLATION_TABLE = "nuclear_translation_table"
    CONFIG_KEY_VCF = "ensembl_vcf_proteindb"
    MITO_TRANSLATION_TABLE = "mito_translation_table"
    HEADER_VAR_PREFIX = "var_prefix"
    REPORT_REFERENCE_SEQ = "report_ref_seq"
    PROTEIN_DB_OUTPUT = "protendb_output_file"
    ANNOTATION_FIELD_NAME = "annotation_field_name"
    AF_FIELD = "af_field"
    AF_THRESHOLD = "af_threshold"
    TRANSCRIPT_INDEX = "transcript_index"
    CONSEQUENCE_INDEX = "consequence_index"
    EXCLUDE_BIOTYPES = "exclude_biotypes"
    EXCLUDE_CONSEQUENCES = "exclude_consequences"
    SKIP_INCLUDING_ALL_CDS = "skip_including_all_cds"
    INCLUDE_BIOTYPES = "include_biotypes"
    INCLUDE_CONSEQUENCES = "include_consequences"
    BIOTYPE_STR = "biotype_str"

    CONFIG_KEY_DATA = "enmsembl_translation"
    CONFIG_TRANSLATION_TABLE = "translation_table"

    def __init__(self, config_file, pipeline_arguments):
        """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
        super(EnsemblDataService, self).__init__(self.CONFIG_KEY_DATA, config_file,
                                                 pipeline_arguments)

        self._proteindb_output = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.PROTEIN_DB_OUTPUT]
        if self.PROTEIN_DB_OUTPUT in self.get_pipeline_parameters():
            self._proteindb_output = self.get_pipeline_parameters()[self.PROTEIN_DB_OUTPUT]

        # Todo: Review if the translation table use in the 3frame method is the same that mito or nuclear?
        self._translation_table = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_TRANSLATION_TABLE]
        if self.CONFIG_TRANSLATION_TABLE in self.get_pipeline_parameters():
            self._translation_table = self.get_pipeline_parameters()[self.CONFIG_TRANSLATION_TABLE]

        self._nuclear_translation_table = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.NUCLEAR_TRANSLATION_TABLE]
        if self.NUCLEAR_TRANSLATION_TABLE in self.get_pipeline_parameters():
            self._nuclear_translation_table = self.get_pipeline_parameters()[self.NUCLEAR_TRANSLATION_TABLE]

        self._mito_translation_table = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.MITO_TRANSLATION_TABLE]
        if self.MITO_TRANSLATION_TABLE in self.get_pipeline_parameters():
            self._mito_translation_table = self.get_pipeline_parameters()[self.MITO_TRANSLATION_TABLE]

        self._header_var_prefix = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.HEADER_VAR_PREFIX]
        if self.HEADER_VAR_PREFIX in self.get_pipeline_parameters():
            self._header_var_prefix = self.get_pipeline_parameters()[self.HEADER_VAR_PREFIX]

        self._report_reference_seq = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.REPORT_REFERENCE_SEQ]
        if self.REPORT_REFERENCE_SEQ in self.get_pipeline_parameters():
            self._report_reference_seq = self.get_pipeline_parameters()[self.REPORT_REFERENCE_SEQ]

        self._annotation_field_name = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.ANNOTATION_FIELD_NAME]
        if self.ANNOTATION_FIELD_NAME in self.get_pipeline_parameters():
            self._annotation_field_name = self.get_pipeline_parameters()[self.ANNOTATION_FIELD_NAME]

        self._af_field = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.AF_FIELD]
        if self.AF_FIELD in self.get_pipeline_parameters():
            self._af_field = self.get_pipeline_parameters()[self.AF_FIELD]

        self._af_threshold = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.AF_THRESHOLD]
        if self.AF_THRESHOLD in self.get_pipeline_parameters():
            self._af_threshold = self.get_pipeline_parameters()[self.AF_THRESHOLD]

        self._transcript_index = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.TRANSCRIPT_INDEX]
        if self.TRANSCRIPT_INDEX in self.get_pipeline_parameters():
            self._transcript_index = self.get_pipeline_parameters()[self.TRANSCRIPT_INDEX]

        self._consequence_index = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.CONSEQUENCE_INDEX]
        if self.CONSEQUENCE_INDEX in self.get_pipeline_parameters():
            self._consequence_index = self.get_pipeline_parameters()[self.CONSEQUENCE_INDEX]

        self._exclude_biotypes = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.EXCLUDE_BIOTYPES]
        if self.EXCLUDE_BIOTYPES in self.get_pipeline_parameters():
            self._exclude_biotypes = self.get_pipeline_parameters()[self.EXCLUDE_BIOTYPES]

        self._exclude_consequences = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.EXCLUDE_CONSEQUENCES]
        if self.EXCLUDE_CONSEQUENCES in self.get_pipeline_parameters():
            self._exclude_consequences = self.get_pipeline_parameters()[self.EXCLUDE_CONSEQUENCES]

        self._skip_including_all_cds = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.SKIP_INCLUDING_ALL_CDS]
        if self.ANNOTATION_FIELD_NAME in self.get_pipeline_parameters():
            self._annotation_field_name = self.get_pipeline_parameters()[self.ANNOTATION_FIELD_NAME]

        self._include_biotypes = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.INCLUDE_BIOTYPES]
        if self.INCLUDE_BIOTYPES in self.get_pipeline_parameters():
            self._include_biotypes = self.get_pipeline_parameters()[self.INCLUDE_BIOTYPES]

        self._include_consequences = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.INCLUDE_CONSEQUENCES]
        if self.INCLUDE_CONSEQUENCES in self.get_pipeline_parameters():
            self._include_consequences = self.get_pipeline_parameters()[self.INCLUDE_CONSEQUENCES]

        self._biotype_str = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
            self.BIOTYPE_STR]
        if self.BIOTYPE_STR in self.get_pipeline_parameters():
            self._biotype_str = self.get_pipeline_parameters()[self.BIOTYPE_STR]

        # Check if some of the variables are pass by commandline

    def three_frame_translation(self, input_file):
        """
        This function translate a transcriptome into a 3'frame translation protein sequence database
        :param input_file: input file
        :param output_file: output file
        :return:
        """

        input_handle = open(input_file, 'r')
        output_handle = open(self._proteindb_output, 'rw')

        for record in SeqIO.parse(input_handle, 'fasta'):
            seq = record.seq
            RF1 = seq.translate(table=self._translation_table)
            RF2 = seq[1::].translate(table=self._translation_table)
            RF3 = seq[2::].translate(table=self._translation_table)

            if record.id == "":
                print("skip entries without id", record.description)
                continue
            output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF1', RF1))
            output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF2', RF2))
            output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF3', RF3))

    @staticmethod
    def check_overlap(var_start, var_end, features_info=[[0, 1, 'type']]):
        """
        This function returns true when the variant overlaps any of the features
        :param var_start: Start location
        :param var_end: End location
        :param features_info: Feature information
        :return:
        """
        if var_start == -1:
            return True
        # check if the var overlaps any of the features
        for feature_pos in features_info:
            pep_start = feature_pos[0]
            pep_end = feature_pos[1]
            if var_start <= pep_start <= var_end:  # fully contained or partial overlap from the end
                return True
            elif var_start <= pep_end <= var_end:  # partial overlap in the begining
                return True
            elif pep_start <= var_start and pep_end >= var_end:  # fully covered
                return True
        return False

    @staticmethod
    def generate_fasta_from_gtf(gtf_file, genome_fasta, fasta_output):
        """use gffread to convert transcript coordinates to genome fasta sequences"""
        if os.path.isfile(fasta_output):
            return fasta_output

        # Todo: @Husen, we should remove calls to external tools.
        cmd = "gffread -F -w {fasta_out} -g {genome_fasta} {gtf_file}".format(
            fasta_out=fasta_output, genome_fasta=genome_fasta, gtf_file=gtf_file)
        os.system(cmd)
        return fasta_output

    @staticmethod
    def get_altseq(ref_seq, ref_allele, var_allele, var_pos, strand, features_info, cds_info=[]):
        """
        the given sequence in the fasta file represents all exons of the transcript combined.
        for protein coding genes, the CDS is specified therefore the sequence position has to be
        calculated based on the CDS's positions
        However, for non-protein coding genes, the whole sequence is used
        """
        alt_seq = ""
        if len(cds_info) == 2:
            start_coding_index = cds_info[0] - 1  # it should be index not pos
            stop_coding_index = cds_info[1]  # get end position of the  last cds
        else:
            start_coding_index = 0
            total_len = 0
            for x in features_info:
                total_len += x[1] - x[0] + 1
            stop_coding_index = total_len  # the features are sorted by end therefroe the end pos of the last item is the last coding nc

        if strand == '-':  # ge the correct orientation, because exons are oredered based on their position
            ref_seq = ref_seq[
                      ::-1]  # in order to calculate from the first base of the first feature (sorted by genomic coordinates)
            ref_allele = ref_allele.complement()  # the reverse will be done on return
            var_allele = var_allele.complement()  # the reverse will be done on return

        ref_seq = ref_seq[
                  start_coding_index:stop_coding_index]  # just keep the coding regions (mostly effective in case of protein-coding genes)
        nc_index = 0
        # Todo: The feature_len variable is not use in the method.
        feature_len = 0
        if len(ref_allele) == len(var_allele) or ref_allele[0] == var_allele[0]:
            for feature in features_info:  # for every exon, cds or stop codon
                if var_pos in range(feature[0], feature[
                                                    1] + 1):  # get index of the var relative to the position of the overlapping feature in the coding region
                    var_index_in_cds = nc_index + (var_pos - feature[0])
                    # modify the coding reference sequence accoding to the var_allele
                    c = len(ref_allele)
                    alt_seq = ref_seq[0:var_index_in_cds] + var_allele + ref_seq[
                                                                         var_index_in_cds + c::]  # variant and ref strand??
                    if strand == '-':
                        return ref_seq[::-1], alt_seq[::-1]
                    else:
                        return ref_seq, alt_seq

                feature_len = (feature[1] - feature[0] + 1)
                nc_index += feature_len
                # feature_seq = coding_ref_seq[nc_index:nc_index+feature_len]

        return ref_seq, alt_seq

    @staticmethod
    def parse_gtf(gtf_fn, gtf_db_file):
        """
        Convert GTF file into a FeatureDB
        :param gtf_fn:
        :param gtf_db_file:
        :return:
        """
        try:
            gffutils.create_db(gtf_fn, gtf_db_file, merge_strategy="create_unique",
                               keep_order=True, disable_infer_transcripts=True, disable_infer_genes=True,
                               verbose=True,
                               force=False)
        except:  # already exists
            pass

        db = gffutils.FeatureDB(gtf_db_file)
        return db

    @staticmethod
    def get_features(db, feature_id, biotype_str, feature_types=['exon']):
        """
        Get chr, genomic positions, strand and biotype for feature_id
        also genomic positions for all its elements (exons/cds&start_codon)
        :param db:
        :param feature_id:
        :param biotype_str:
        :param feature_types:
        :return:
        """
        feature = db[feature_id]
        coding_features = []
        features = db.children(feature, featuretype=feature_types, order_by='end')
        for f in features:
            f_type = f.featuretype
            coding_features.append([f.start, f.end, f_type])
        return feature.chrom, feature.strand, coding_features, feature.attributes[biotype_str][0]

    @staticmethod
    def get_orfs(ref_seq, alt_seq, trans_table, num_orfs=1):
        """
        Translate the coding_ref and the coding_alt into ORFs
        :param ref_seq:
        :param alt_seq:
        :param trans_table:
        :param num_orfs:
        :return:
        """

        ref_orfs = []
        alt_orfs = []
        for n in range(0, num_orfs):
            ref_orfs.append(ref_seq[n::].translate(trans_table))
            alt_orfs.append(alt_seq[n::].translate(trans_table))

        return ref_orfs, alt_orfs

    def vcf_to_proteindb(self, vcf_file, genome_fasta, gtf_fasta):
        """
        Generate peps for variants by modifying sequences of affected transcripts (VEP annotated).
        It only considers variants within potential coding regions of the transcript
        (CDSs & stop codons for protein-coding genes, exons for non-protein coding genes).
        :param vcf_file:
        :param genome_fasta:
        :param gtf_fasta:
        :return:
        """
        transcripts_dict = SeqIO.index(genome_fasta, "fasta")
        # handle cases where the transript has version in the GTF but not in the VCF
        transcript_id_mapping = {k.split('.')[0]: k for k in transcripts_dict.keys()}

        with open(self._proteindb_output, 'w') as prots_fn:
            vcf_reader = vcf.Reader(open(vcf_file, 'r'))
            for record in vcf_reader:
                # only process variants above a given allele frequencey threshold
                # get AF from the INFO field
                try:
                    af = float(record.INFO[self._af_field])
                except TypeError:
                    af = float(record.INFO[self._af_field][0])
                except KeyError:
                    continue

                # check if the AF passed the threshold
                if af < self._af_threshold:
                    continue

                trans_table = self._nuclear_translation_table
                consequences = []
                if str(record.CHROM).lstrip('chr').upper() in ['M', 'MT']:
                    trans_table = self._mito_translation_table

                processed_transcript_allele = []
                for transcript_record in record.INFO[self._annotation_field_name]:
                    transcript_info = transcript_record.split('|')

                    try:
                        consequence = transcript_info[self._consequence_index]
                    except IndexError:
                        print("Give a valid index for the consequence in the INFO field for: ", transcript_record)
                        continue
                    consequences.append(consequence)

                    try:
                        transcript_id = transcript_info[self._transcript_index]
                    except IndexError:
                        print("Give a valid index for the Transcript ID in the INFO field for: ", transcript_record)
                        continue
                    if transcript_id == "":
                        continue

                    try:
                        transcript_id_v = transcript_id_mapping[transcript_id]
                    except KeyError:
                        transcript_id_v = transcript_id

                    try:
                        row = transcripts_dict[transcript_id_v]
                        ref_seq = row.seq  # get the seq and desc for the transcript from the fasta of the gtf
                        desc = str(row.description)
                    except KeyError:
                        print("Transcript {} not found in fasta of the GTF file {}".format(transcript_id_v, record))
                        continue

                    feature_types = ['exon']
                    # check if cds info exists in the fasta header otherwise translate all exons
                    cds_info = []
                    num_orfs = 3
                    if 'CDS=' in desc:
                        try:
                            cds_info = [int(x) for x in desc.split(' ')[1].split('=')[1].split('-')]
                            feature_types = ['CDS', 'stop_codon']
                            num_orfs = 1
                        except (ValueError, IndexError):
                            print("Could not extra cds position from fasta header for: ", desc)
                            pass

                    chrom, strand, features_info, feature_biotype = get_features(gtf_fasta, transcript_id_v, self._biotype_str, feature_types)

                    # skip transcripts with unwanted consequences
                    if (consequence in self._exclude_consequences or
                            (consequence not in self._include_consequences and self._include_consequences != ['all'])
                    ):
                        continue

                    # only include features that have the specified biotypes or they have CDSs info
                    if 'CDS' in feature_types and not self._skip_including_all_cds:
                        pass
                    elif (feature_biotype in self._exclude_biotypes or
                          (feature_biotype not in self._include_biotypes and self._include_biotypes != ['all'])):
                        continue

                    for alt in record.ALT:  # in cases of multiple alternative alleles consider all
                        if transcript_id + str(
                                alt) not in processed_transcript_allele:  # because VEP reports affected transcripts per alt allele
                            processed_transcript_allele.append(transcript_id + str(alt))
                            "for non-CDSs, only consider the exon that actually overlaps the variant"
                            if (chrom.lstrip("chr") == str(record.CHROM).lstrip("chr") and
                                    check_overlap(record.POS, record.POS + len(alt), features_info)):
                                coding_ref_seq, coding_alt_seq = get_altseq(ref_seq, Seq(str(record.REF)),
                                                                            Seq(str(alt)),
                                                                            int(record.POS), strand, features_info,
                                                                            cds_info)
                                if coding_alt_seq != "":
                                    ref_orfs, alt_orfs = get_orfs(coding_ref_seq, coding_alt_seq, trans_table, num_orfs)

                                    record_id = ""
                                    if record.ID:
                                        record_id = '_' + str(record.ID)
                                    self.write_output(seq_id='_'.join([self._header_var_prefix + str(record_id), '.'.join([str(record.CHROM), str(record.POS), str(record.REF), str(alt)]),
                                                                  transcript_id_v]),
                                                 desc=feature_biotype + ":" + consequence,
                                                 seqs=alt_orfs,
                                                 prots_fn=prots_fn)

                                    if self._report_reference_seq:
                                        self.write_output(seq_id=transcript_id_v,
                                                     desc=feature_biotype,
                                                     seqs=ref_orfs,
                                                     prots_fn=prots_fn)

        return self._proteindb_output

    def write_output(self, seq_id, desc, seqs, prots_fn):
        "write the orfs to the output file"
        for i, alt_orf in enumerate(seqs):
            orf_num = ""
            if i > 0:  # only add _num when multiple ORFs are generated (e.g in 3 ORF)
                orf_num = "_" + str(i + 1)

            prots_fn.write('>{} {}\n{}\n'.format(seq_id + orf_num, desc, alt_orf))


if __name__ == '__main__':
    print("ERROR: This script is part of a pipeline collection and it is not meant to be run in stand alone mode")
