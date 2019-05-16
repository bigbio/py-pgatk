import re

from Bio import SeqIO

from pypgatk.cgenomes.models import SNP
from pypgatk.toolbox.general import ParameterConfiguration


class CancerGenomesService(ParameterConfiguration):
    CONFIG_CANCER_GENOMES_MUTATION_FILE = 'mutation_file'
    CONFIG_COMPLETE_GENES_FILE = "genes_file"
    CONFIG_OUTPUT_FILE = "output_file"
    CONFIG_COSMIC_DATA = "cosmic_data"

    def __init__(self, config_file, pipeline_arguments):
        """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
        super(CancerGenomesService, self).__init__(self.CONFIG_COSMIC_DATA, config_file, pipeline_arguments)

        if self.CONFIG_CANCER_GENOMES_MUTATION_FILE in self.get_pipeline_parameters():
            self._local_mutation_file = self.get_pipeline_parameters()[self.CONFIG_CANCER_GENOMES_MUTATION_FILE]

        if self.CONFIG_COMPLETE_GENES_FILE in self.get_pipeline_parameters():
            self._local_complete_genes = self.get_pipeline_parameters()[self.CONFIG_COMPLETE_GENES_FILE]

        if self.CONFIG_OUTPUT_FILE in self.get_pipeline_parameters():
            self._local_output_file = self.get_pipeline_parameters()[self.CONFIG_OUTPUT_FILE]

    def cosmic_to_proteindb(self):
        """
        This function translate the mutation file + COSMIC genes into a protein Fasta database. The
        method write into the file system the output Fasta.
        :return:
        """
        self.get_logger().debug("Starting reading the All cosmic genes")
        COSMIC_CDS_DB = SeqIO.index(self._local_complete_genes, 'fasta')  # All_COSMIC_Genes.fasta

        cosmic_input = open(self._local_mutation_file, 'r')  # CosmicMutantExport.tsv

        header = cosmic_input.readline().split("\t")

        gene_col = header.index("Gene name")
        enst_col = header.index("Accession Number")
        cds_col = header.index("Mutation CDS")
        aa_col = header.index("Mutation AA")
        muttype_col = header.index("Mutation Description")

        output = open(self._local_output_file, 'w')

        mutation_dic = {}
        nucleotide = ["A", "T", "C", "G"]
        self.get_logger().debug("Reading input CosmicMutantExport.tsv ...")
        line_counter = 1
        for line in cosmic_input:
            if line_counter % 10000 == 0:
                msg = "Number of lines finished -- '{}'".format(line_counter)
                self.get_logger().debug(msg)
            line_counter += 1
            row = line.strip().split("\t")
            if "coding silent" in row[muttype_col]:
                continue

            snp = SNP(gene=row[gene_col], mrna=row[enst_col], dna_mut=row[cds_col], aa_mut=row[aa_col],
                      type=row[muttype_col])
            header = "COSMIC:%s:%s:%s" % (snp.gene, snp.aa_mut, snp.type.replace(" ", ""))
            try:
                seq = COSMIC_CDS_DB[snp.gene].seq
            except KeyError:  # geneID is not in All_COSMIC_Genes.fasta
                continue

            mut_pro_seq = ""
            if "?" not in row[cds_col]:  # unambiguous DNA change known in CDS sequence
                positions = re.findall(r'\d+', snp.dna_mut)
                if ">" in snp.dna_mut and len(positions) == 1:  # Substitution
                    tmplist = snp.dna_mut.split(">")
                    ref_dna = re.sub("[^A-Z]+", "", tmplist[0])
                    mut_dna = re.sub("[^A-Z]+", "", tmplist[1])
                    index = int(positions[0]) - 1

                    if ref_dna == str(seq[index]).upper() and mut_dna in nucleotide:  #
                        seq_mut = seq[:index] + mut_dna + seq[index + 1:]
                        mut_pro_seq = seq_mut.translate(to_stop=True)
                elif "ins" in snp.dna_mut:
                    index = snp.dna_mut.index("ins")
                    insert_dna = snp.dna_mut[index + 3:]
                    if insert_dna.isalpha():
                        print(insert_dna, snp.dna_mut, snp.aa_mut)
                        ins_index1 = int(positions[0])
                        seq_mut = seq[:ins_index1] + insert_dna + seq[ins_index1:]
                        mut_pro_seq = seq_mut.translate(to_stop=True)

                elif "del" in snp.dna_mut:
                    if len(positions) == 2:
                        del_index1 = int(positions[0]) - 1
                        del_index2 = int(positions[1])
                        seq_mut = seq[:del_index1] + seq[del_index2:]
                        mut_pro_seq = seq_mut.translate(to_stop=True)
                    elif len(positions) == 1:
                        del_index1 = int(positions[0]) - 1
                        seq_mut = seq[:del_index1] + seq[del_index1 + 1:]
                        mut_pro_seq = seq_mut.translate(to_stop=True)
            else:
                if "?" not in row[aa_col]:  # unambiguous aa change known in protein sequence
                    positions = re.findall(r'\d+', snp.aa_mut)
                    protein_seq = str(seq.translate(to_stop=True))

                    if "Missense" in snp.type:
                        mut_aa = snp.aa_mut[-1]
                        index = int(positions[0]) - 1
                        mut_pro_seq = protein_seq[:index] + mut_aa + protein_seq[index + 1:]
                    elif "Nonsense" in snp.type:
                        index = int(positions[0]) - 1
                        mut_pro_seq = protein_seq[:index]
                    elif "Insertion - In frame" in snp.type:
                        index = snp.aa_mut.index("ins")
                        insert_aa = snp.aa_mut[index + 3:]
                        if insert_aa.isalpha():
                            ins_index1 = int(positions[0])
                            mut_pro_seq = protein_seq[:ins_index1] + insert_aa + protein_seq[ins_index1:]
                    elif "Deletion - In frame" in snp.type:
                        if len(positions) == 2:
                            del_index1 = int(positions[0]) - 1
                            del_index2 = int(positions[1])
                            mut_pro_seq = protein_seq[:del_index1] + protein_seq[del_index2:]
                        elif len(positions) == 1:
                            del_index1 = int(positions[0]) - 1
                            mut_pro_seq = protein_seq[:del_index1] + protein_seq[del_index1 + 1:]
                    elif "Complex" in snp.type and "frameshift" not in snp.type:
                        try:
                            index = snp.aa_mut.index(">")
                        except ValueError:
                            # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                            continue
                        mut_aa = snp.aa_mut[index + 1:]
                        if "deletion" in snp.type:
                            try:
                                del_index1 = int(positions[0]) - 1
                                del_index2 = int(positions[1])
                                mut_pro_seq = protein_seq[:del_index1] + mut_aa + protein_seq[del_index2:]
                            except IndexError:
                                # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                                continue
                        elif "insertion" in snp.type:
                            try:
                                ins_index1 = int(positions[0]) - 1
                            except IndexError:
                                # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                                continue
                            mut_pro_seq = protein_seq[:ins_index1] + mut_aa + protein_seq[ins_index1 + 1:]
                        elif "compound substitution" in snp.type:
                            if "*" not in mut_aa:
                                try:
                                    del_index1 = int(positions[0]) - 1
                                    del_index2 = int(positions[1])
                                    mut_pro_seq = protein_seq[:del_index1] + mut_aa + protein_seq[del_index2:]
                                except IndexError:
                                    # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                                    continue
                            else:
                                try:
                                    del_index1 = int(positions[0]) - 1
                                    mut_pro_seq = protein_seq[:del_index1] + mut_aa.replace("*", "")
                                except IndexError:
                                    # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                                    continue
            if mut_pro_seq != "":
                entry = ">%s\n%s\n" % (header, mut_pro_seq)
                if header not in mutation_dic:
                    output.write(entry)
                    mutation_dic[header] = 1

        self.get_logger().debug("COSMIC contains in total", len(mutation_dic), "non redundant mutations")
        cosmic_input.close()
        output.close()

    def cbioportal_to_proteindb(self):

        mutfile = open(self._local_mutation_file, "r")
        fafile = SeqIO.parse(self._local_complete_genes, "fasta")
        output = open(self._local_output_file, "w")

        seq_dic = {}
        for record in fafile:
            newacc = record.id.split(".")[0]
            if newacc not in seq_dic:
                seq_dic[newacc] = record.seq

        print(len(seq_dic))

        f1 = mutfile.readline()
        if f1[0] == "#":
            header = mutfile.readline().strip().split("\t")
        else:
            header = f1.strip().split("\t")

        # print(header)
        pos_col = header.index("HGVSc")
        enst_col = header.index("Transcript_ID")

        class_col = header.index("Variant_Classification")
        type_col = header.index("Variant_Type")
        aa_col = header.index("HGVSp_Short")

        nucleotide = ["A", "T", "C", "G"]
        mutclass = ["Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                    "Nonsense_Mutation"]

        for line in mutfile:
            row = line.strip().split("\t")
            gene = row[0]
            try:
                pos = row[pos_col]
                enst = row[enst_col]

                seq_mut = ""
                aa_mut = row[aa_col]

                vartype = row[type_col]
                varclass = row[class_col]
            except IndexError:
                print(row)
                continue

            if varclass not in mutclass:
                continue

            if enst in seq_dic:
                seq = seq_dic[enst]
            else:
                print("%s not found" % enst)
                continue

            if ":" in pos:
                cdna_pos = pos.split(":")[1]
            else:
                cdna_pos = pos

            if vartype == "SNP":
                enst_pos = int(re.findall(r'\d+', cdna_pos)[0])
                idx = pos.index(">")
                ref_dna = pos[idx - 1]
                mut_dna = pos[idx + 1]

                if mut_dna not in nucleotide:
                    print(mut_dna, "is not a nucleotide base", pos)
                    continue
                try:
                    if ref_dna == seq[enst_pos - 1]:
                        seq_mut = seq[:enst_pos - 1] + mut_dna + seq[enst_pos:]
                    else:
                        print("incorrect substitution, unmatched nucleotide", pos, enst)
                except IndexError:
                    print("incorrect substitution, out of index", pos)
            elif vartype == "DEL":
                try:
                    enst_pos = int(re.findall(r'\d+', cdna_pos.split("_")[0])[0])
                except IndexError:
                    print("incorrect del format", pos)
                    continue
                del_dna = pos.split("del")[1]
                if del_dna == seq[enst_pos - 1:enst_pos - 1 + len(del_dna)]:
                    seq_mut = seq[:enst_pos - 1] + seq[enst_pos - 1 + len(del_dna):]
                else:
                    print("incorrect deletion, unmatched nucleotide", pos)

            elif vartype == "INS":
                enst_pos = int(re.findall(r'\d+', cdna_pos.split("_")[0])[0])
                if "ins" in pos:
                    ins_dna = pos.split("ins")[1]
                elif "dup" in pos:
                    ins_dna = pos.split("dup")[1]
                    if len(ins_dna) > 1:
                        enst_pos = int(re.findall(r'\d+', cdna_pos.split("_")[1])[0])
                else:
                    print("unexpected insertion format")
                    continue;

                seq_mut = seq[:enst_pos] + ins_dna + seq[enst_pos:]

            if seq_mut == "":
                continue;

            mut_pro_seq = seq_mut.translate(to_stop=True)
            if len(mut_pro_seq) > 6:
                header = "Mutation:%s:%s:%s:%s" % (enst, gene, aa_mut, varclass)
                output.write(">%s\n%s\n" % (header, mut_pro_seq))

        output.close()
        mutfile.close()
        fafile.close()
