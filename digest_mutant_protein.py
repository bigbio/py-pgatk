import sys
import re
from Bio import SeqIO
from collections import OrderedDict

def TRYPSIN(proseq,miss_cleavage):
    peptides=[]
    cut_sites=[0]
    for i in range(0,len(proseq)-1):
        if proseq[i]=='K' and proseq[i+1]!='P':
            cut_sites.append(i+1)
        elif proseq[i]=='R' and proseq[i+1]!='P':
            cut_sites.append(i+1)
    
    if cut_sites[-1]!=len(proseq):
        cut_sites.append(len(proseq))
    
    if  miss_cleavage==0:
        for j in range(0,len(cut_sites)-1):
            peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
    
    elif miss_cleavage==1:
        for j in range(0,len(cut_sites)-2):
            peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
            peptides.append(proseq[cut_sites[j]:cut_sites[j+2]])
            
        peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])
    
    elif miss_cleavage==2:
        for j in range(0,len(cut_sites)-3):
            peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
            peptides.append(proseq[cut_sites[j]:cut_sites[j+2]])
            peptides.append(proseq[cut_sites[j]:cut_sites[j+3]])
        
        peptides.append(proseq[cut_sites[-3]:cut_sites[-2]])
        peptides.append(proseq[cut_sites[-3]:cut_sites[-1]])
        peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])
    
    return peptides

handle1=SeqIO.parse(sys.argv[1],'fasta') # All_COSMIC_Genes.fasta


filelist=sys.argv[3:]
handle_list=[]
for f in filelist:
    handle_list.append(SeqIO.parse(f,'fasta')) # Cosmic_mutant_proteins.fasta

output=open(sys.argv[2],'w')

peptidome={}

for record in handle1:
    cds_seq = record.seq
    aa_seq = cds_seq.translate(to_stop=True,stop_symbol="")
    peptide_list=TRYPSIN(str(aa_seq),0)
    for peptide in peptide_list:
        if len(peptide) in range(6,41):
            if peptide not in peptidome:
                peptidome[peptide.replace("I","L")]=1

print("known peptides number",len(peptidome))
handle1.close()

var_peptidome=OrderedDict()

for h in handle_list:
    for record in h:
        proseq=record.seq
        descrip=record.description
        if len(proseq)>5:
            peptide_list=TRYPSIN(str(proseq),0)
            for peptide in peptide_list:
                if len(peptide) in range(6,41):
                    peptide1=peptide.replace("I","L")
                    if peptide1 not in peptidome:
                        des_list=descrip.split(":")
                        des_list[0]="COSMIC"
                        type=des_list[-1]
                        snp=des_list[-2]
                        if "Missense" in type:
                            try:
                                mut_pos=int(re.findall(r'\d+',snp)[0])
                                index=str(proseq).index(peptide)
                                pos=mut_pos-index
                                des_list.append(str(pos))
                            except IndexError:
                                continue;
                    
                        new_description=":".join(des_list).replace("*","-")
                        if peptide not in var_peptidome:
                            var_peptidome[peptide]=[new_description]
                        else:
                            var_peptidome[peptide].append(new_description)
    h.close()
    print("file process done")
    

print("mut peptide numbers",len(var_peptidome))
for pep in var_peptidome.keys():
    acc=";".join(set(var_peptidome[pep]))
    output.write(">%s\n%s\n" % (acc,pep))

output.close()
