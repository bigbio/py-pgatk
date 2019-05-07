import sys
from Bio import Seq,SeqIO
from Bio.Alphabet import IUPAC
import os
import getopt

table = 1

if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command, please read the mannual")
    print("Example: python 3frame_translation.py --input transcriptome.fasta --output transcriptome.3FT.fasta --translate_table 1")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=',
                                                         'output=',
                                                         'translate_table='])                                                
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--output': output_file=arg
        elif opt == '--translate_table': table=int(arg)
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()


print("translation table %d is used" % table)
print("find the correct code at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")

handle=open(input_file,'r')
output=open(output_file,'w')

for record in SeqIO.parse(handle,'fasta'):
    seq=record.seq
    RF1=seq.translate(table = table)
    RF2=seq[1::].translate(table = table)
    RF3=seq[2::].translate(table = table)
    
    if record.id=="":
        print("skip entries without id",record.description)
        continue
    output.write("%s\n%s\n" % ('>'+record.id+'_RF1',RF1))
    output.write("%s\n%s\n" % ('>'+record.id+'_RF2',RF2))
    output.write("%s\n%s\n" % ('>'+record.id+'_RF3',RF3))
            
handle.close()
output.close()
