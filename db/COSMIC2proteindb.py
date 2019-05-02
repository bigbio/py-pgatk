'''
    this script is make a fasta file of mutant protein sequence from COSMIC database.
    It requires two input files: CosmicMutantExport.tsv and All_COSMIC_Genes.fasta from http://cancer.sanger.ac.uk/cosmic/download
'''

import getopt
import re
import sys

from Bio import SeqIO

if len(sys.argv[1:]) <= 1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command!")
    print(
        "Example: python3 COSMIC2proteindb.py --input CosmicMutantExport.tsv --fa All_COSMIC_Genes.fasta --output mutproteins.fa")
else:
    options, remainder = getopt.getopt(sys.argv[1:], '', ['input=', 'fa=', 'output='])
    for opt, arg in options:
        if opt == '--input':
            input_file = arg
        elif opt == '--output':
            output_file = arg
        elif opt == '--fa':
            fa_file = arg
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt);
            sys.exit()


class SNP(object):
    def __init__(self, gene=None, mRNA=None, dna_mut=None, aa_mut=None, type=None):
        self.gene = gene
        self.mRNA = mRNA
        self.aa_mut = aa_mut
        self.type = type
        self.dna_mut = dna_mut


print("reading gene CDS sequences ...")
COSMIC_CDS_DB = SeqIO.index(fa_file, 'fasta')  # All_COSMIC_Genes.fasta

cosmic_input = open(input_file, 'r')  # CosmicMutantExport.tsv

header = cosmic_input.readline().split("\t")

gene_col = header.index("Gene name")
ENST_col = header.index("Accession Number")
cds_col = header.index("Mutation CDS")
aa_col = header.index("Mutation AA")
muttype_col = header.index("Mutation Description")

output = open(output_file, 'w')

mutation_dic = {}
nucleotide = ["A", "T", "C", "G"]
print("reading input CosmicMutantExport.tsv ...")
line_counter = 1
for line in cosmic_input:
    if line_counter % 10000 == 0:
        print("Number of lines finished", line_counter)
    line_counter += 1
    row = line.strip().split("\t")
    if "coding silent" in row[muttype_col]:
        continue;

    snp = SNP(gene=row[gene_col], mRNA=row[ENST_col], dna_mut=row[cds_col], aa_mut=row[aa_col], type=row[muttype_col])
    header = "COSMIC:%s:%s:%s" % (snp.gene, snp.aa_mut, snp.type.replace(" ", ""))
    try:
        seq = COSMIC_CDS_DB[snp.gene].seq
    except KeyError:  # geneID is not in All_COSMIC_Genes.fasta
        continue;

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
            index = snp.dna_mut.index("del")
            if len(positions) > 1:
                del_index1 = int(positions[0]) - 1
                del_index2 = int(positions[1])
                seq_mut = seq[:del_index1] + seq[del_index2:]
                mut_pro_seq = seq_mut.translate(to_stop=True)
            else:
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
                try:
                    index = snp.aa_mut.index("del")
                except ValueError:
                    # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type
                    continue;
                if len(positions) > 1:
                    del_index1 = int(positions[0]) - 1
                    del_index2 = int(positions[1])
                    mut_pro_seq = protein_seq[:del_index1] + protein_seq[del_index2:]
                else:
                    del_index1 = int(positions[0]) - 1
                    mut_pro_seq = protein_seq[:del_index1] + protein_seq[del_index1 + 1:]
            elif "Complex" in snp.type and "frameshift" not in snp.type:
                try:
                    index = snp.aa_mut.index(">")
                except ValueError:
                    # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                    continue;
                mut_aa = snp.aa_mut[index + 1:]
                if "deletion" in snp.type:
                    try:
                        del_index1 = int(positions[0]) - 1
                        del_index2 = int(positions[1])
                        mut_pro_seq = protein_seq[:del_index1] + mut_aa + protein_seq[del_index2:]
                    except IndexError:
                        # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                        continue;
                elif "insertion" in snp.type:
                    try:
                        ins_index1 = int(positions[0]) - 1
                    except IndexError:
                        # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                        continue;
                    mut_pro_seq = protein_seq[:ins_index1] + mut_aa + protein_seq[ins_index1 + 1:]
                elif "compound substitution" in snp.type:
                    if "*" not in mut_aa:
                        try:
                            del_index1 = int(positions[0]) - 1
                            del_index2 = int(positions[1])
                            mut_pro_seq = protein_seq[:del_index1] + mut_aa + protein_seq[del_index2:]
                        except IndexError:
                            # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                            continue;
                    else:
                        try:
                            del_index1 = int(positions[0]) - 1
                            del_index2 = int(positions[1])
                            mut_pro_seq = protein_seq[:del_index1] + mut_aa.replace("*", "")
                        except IndexError:
                            # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                            continue;
    if mut_pro_seq != "":
        entry = ">%s\n%s\n" % (header, mut_pro_seq)
        if header not in mutation_dic:
            output.write(entry)
            mutation_dic[header] = 1

print("COSMIC contains in total", len(mutation_dic), "non redundant mutations")
cosmic_input.close()
output.close()
