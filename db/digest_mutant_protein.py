import sys
import os
import getopt
import re
from Bio import SeqIO
from collections import OrderedDict

def TRYPSIN(proseq, miss_cleavage):
    peptides = []
    peptide = ''
    for c, aa in enumerate(proseq):
        peptide+=aa
        next_aa = ''
        try:
            next_aa = proseq[c+1]
        except IndexError:
            pass
            
        if aa in ['K', 'R'] and next_aa!= 'P': # for trypsin peptides
            if len(peptide)>0:
                peptides.append(peptide)
            peptide = ''
            continue
            
    if len(peptide)>0:
        peptides.append(peptide)
    
    peptides_with_miss_cleavage = []
    for i in range(1, miss_cleavage+1):
        for j,pep in enumerate(peptides):
            if j+i<len(peptides):
                peptide = ''.join([x for x in (peptides[j:j+i+1])])
                peptides_with_miss_cleavage.append(peptide)
    
    peptides.extend(peptides_with_miss_cleavage)

    return peptides


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command!")
    print("Example: python3 db/mutation2proteindb.py --input mutproteins.1.fa,mutproteins.2.fa --cds Ensembl75+90.human.cds.all.fa --prefix Mutation --output mutpeptides.fa")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=','cds=','output=','prefix='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--output': output_file=arg
        elif opt == '--cds': cds_file=arg
        elif opt == '--prefix': header_prefix=arg
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()

handle1 = SeqIO.parse(cds_file, 'fasta')  # gene cds sequence
filelist = input_file.split(",")
handle_list = []
for f in filelist:
    handle_list.append(SeqIO.parse(f, 'fasta'))  # mutproteins.fasta

output = open(output_file, 'w')

peptidome = {}

for record in handle1:
    cds_seq = record.seq
    aa_seq = cds_seq.translate(to_stop=True, stop_symbol="")
    peptide_list = TRYPSIN(str(aa_seq), 0)
    for peptide in peptide_list:
        if len(peptide) in range(6, 41):
            if peptide not in peptidome:
                peptidome[peptide.replace("I", "L")] = 1

print("known peptides number", len(peptidome))
handle1.close()

var_peptidome = OrderedDict()

for h in handle_list:
    for record in h:
        proseq = record.seq
        descrip = record.description
        if len(proseq) > 5:
            peptide_list = TRYPSIN(str(proseq), 0)
            for peptide in peptide_list:
                if len(peptide) in range(6, 41):
                    peptide1 = peptide.replace("I", "L")
                    if peptide1 not in peptidome:
                        des_list = descrip.split(":")
                        des_list[0] = header_prefix
                        type = des_list[-1]
                        snp = des_list[-2]
                        if "Missense" in type:
                            try:
                                mut_pos = int(re.findall(r'\d+', snp)[0])
                                index = str(proseq).index(peptide)
                                pos = mut_pos - index
                                des_list.append(str(pos))
                            except IndexError:
                                continue;

                        new_description = ":".join(des_list).replace("*", "-")
                        if peptide not in var_peptidome:
                            var_peptidome[peptide] = [new_description]
                        else:
                            var_peptidome[peptide].append(new_description)
    h.close()
    print("file process done")

print("mut peptide numbers", len(var_peptidome))
for pep in var_peptidome.keys():
    acc = ";".join(set(var_peptidome[pep]))
    output.write(">%s\n%s\n" % (acc, pep))

output.close()
