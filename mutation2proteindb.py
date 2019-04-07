import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re

mutfile = open(sys.argv[1],"r")
fafile = SeqIO.parse(sys.argv[2],"fasta")
output = open(sys.argv[3],"w")
    
seqdic={}
for record in fafile:
    newacc=record.id.split(".")[0]
    if newacc not in seqdic:
        seqdic[newacc]=record.seq

print(len(seqdic))

f1 = mutfile.readline()
if f1[0]=="#":
    header = mutfile.readline().strip().split("\t")
else:
    header = f1.strip().split("\t")

print (header)
pos_col = header.index("HGVSc")
enst_col = header.index("Transcript_ID")

class_col = header.index("Variant_Classification")
type_col = header.index("Variant_Type")
aa_col = header.index("HGVSp_Short")


nucleotide=["A","T","C","G"]
mutclass=["Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation"]

for line in mutfile:
    row = line.strip().split("\t")
    gene = row[0]
    try:
        pos = row[pos_col] 
        enst = row[enst_col]
 
        seq_mut=""
        aa_mut = row[aa_col]

        vartype = row[type_col]
        varclass = row[class_col]
    except IndexError:
        print (row)
        continue;

    if varclass not in mutclass:
        continue;
    
    if enst in seqdic:
        seq=seqdic[enst]
    else:
        print("%s not found" % (enst))
        continue;

    if ":" in pos:
        cdna_pos=pos.split(":")[1]
    else:
        cdna_pos=pos
    
    if vartype=="SNP":
        enst_pos= int(re.findall(r'\d+',cdna_pos)[0])
        idx = pos.index(">")
        ref_dna = pos[idx-1]
        mut_dna = pos[idx+1]
       
        if mut_dna not in nucleotide:
            print(mut_dna,"is not a nucleotide base",pos)
            continue;
        try:
            if ref_dna==seq[enst_pos-1]:
                seq_mut=seq[:enst_pos-1]+mut_dna+seq[enst_pos:]   
            else:
                print("incorrect substitution, unmatched nucleotide",pos,enst)
        except IndexError:
            print("incorrect substitution, out of index",pos)
    elif vartype=="DEL":
        try:
            enst_pos = int(re.findall(r'\d+',cdna_pos.split("_")[0])[0])
        except IndexError:
            print("incorrect del format",pos)
            continue;
        del_dna=pos.split("del")[1]
        if del_dna == seq[enst_pos-1:enst_pos-1+len(del_dna)]:
            seq_mut=seq[:enst_pos-1]+seq[enst_pos-1+len(del_dna):]
        else:
            print("incorrect deletion, unmatched nucleotide",pos)
              
    elif vartype=="INS":
        enst_pos= int(re.findall(r'\d+',cdna_pos.split("_")[0])[0])
        if "ins" in pos:
            ins_dna=pos.split("ins")[1]
        elif "dup" in pos:
            ins_dna=pos.split("dup")[1]
            if len(ins_dna)>1:
                enst_pos= int(re.findall(r'\d+',cdna_pos.split("_")[1])[0])
        else:
            print("unexpected insertion format")
            continue;

        seq_mut=seq[:enst_pos]+ins_dna+seq[enst_pos:]           
    
    if seq_mut=="":
        continue;

    mut_pro_seq=seq_mut.translate(to_stop=True)
    if len(mut_pro_seq)>6:
        header="Mutation:%s:%s:%s:%s" % (enst,gene,aa_mut,varclass)
        output.write(">%s\n%s\n" % (header,mut_pro_seq))

output.close()
mutfile.close()
fafile.close()
