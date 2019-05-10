"""
Translate transcripts_fasta
"""

import sys
from Bio import SeqIO
import argparse

def parse_commandline_args():
    
    parser = argparse.ArgumentParser(description='Generate peptides based on DNA variants (VCF)')
    
    #input files
    parser.add_argument('--transcripts_fasta', help='path to the gene annotations file', required=True)
    
    #translation params
    parser.add_argument('--trans_table', default=1, type=int, help='nuclear_trans_table (default 1)')
    parser.add_argument('--num_orfs', default=3, type=int, help='for genes lacking CDS info, specify number of ORFs (default 0)')
    parser.add_argument('--num_orfs_complement', default=0, type=int, help='for genes lacking CDS info, specify number of ORFs from the reverse side (default 0')
    
    #output params
    parser.add_argument('--proteindb_output_file', default = "peps.fa", help="output file name, exits if already exists")
    
    #GTF params and filters
    parser.add_argument('--skip_including_all_CDSs', action='store_const', const=True, help="by default any transcript that has a defined CDS will be used, this option disables this features instead it only depends on the biotypes")
    parser.add_argument('--included_biotypes', nargs='*', default = [], 
                        help="included_biotypes")
    parser.add_argument('--excluded_biotypes', nargs='*', default = [], help="excluded_biotypes")
    parser.add_argument('--biotype_str', default='transcript_biotype', type=str, help='string used to identify gene/transcript biotype in the gtf file.')
    parser.add_argument('--expression_str', default="", type=str, help='string to be used for extracting expression value (TPM, FPKM, etc).')
    parser.add_argument('--expression_thresh', default=5.0, type=float, help='threshold used to filter transcripts based on their expression values')
    
    return parser.parse_args(sys.argv[1:])
    
    
args = parse_commandline_args()
    
def get_orfs(ref_seq, trans_table, num_orfs, num_orfs_complement):
    """translate the coding_ref into ORFs"""
    
    ref_orfs = []
    for n in range(0, num_orfs):
        ref_orfs.append(ref_seq[n::].translate(trans_table))    
    
    rev_ref_seq = ref_seq.reverse_complement()
    for n in range(0, num_orfs_complement):
        ref_orfs.append(rev_ref_seq[n::].translate(trans_table))
    
    return ref_orfs
    
    
def translate_transcript(transcripts_fasta,  
                  proteindb_output_file,  
                  trans_table, num_orfs, num_orfs_complement,
                  biotype_str,
                  excluded_biotypes, included_biotypes,
                  skip_including_all_CDSs,
                  expression_str, expression_thresh
                  ):
    transcripts_dict = SeqIO.index(transcripts_fasta, "fasta")
    
    with open(proteindb_output_file , 'w') as prots_fn:
        for transcript_id in transcripts_dict.keys():
            
            ref_seq = transcripts_dict[transcript_id].seq#get the seq and desc for the transcript from the fasta of the gtf
            desc = str(transcripts_dict[transcript_id].description)
            
            key_values = {} #extract key=value in the desc into a dict
            for value in desc.split(' '):
                try:
                    key_values[value.split('=')[0]]=value.split('=')[1]
                except IndexError:
                    continue
            this_num_orfs = num_orfs
            this_num_orfs_complement = num_orfs_complement
            
            feature_biotype = ""
            try:
                feature_biotype = key_values[biotype_str]
            except KeyError:
                print("Biotype info was not found in the header using {} for record {} {}".format(
                    biotype_str, transcript_id, desc))
            
            #only include features that have the specified biotypes or they have CDSs info
            if 'CDS' in key_values.keys() and not skip_including_all_CDSs:
                pass
            elif feature_biotype=="" or (feature_biotype in excluded_biotypes or 
                (feature_biotype not in included_biotypes and included_biotypes!=['all'])):
                continue
            
            #check wether to filter on expression and if it passes
            if expression_str:
                try:
                    if float(key_values[expression_str]) < expression_thresh:
                        continue
                except KeyError:
                    print("Expression information not found in the fasta header with expression_str: {} for record {} {}".format(expression_str, transcript_id, desc))
                    continue
                except TypeError:
                    print("Expression value is not of valid type (float) at record: {} {}".format(transcript_id, key_values[expression_str]))
                    continue
            
            #check if cds info exists in the fasta header otherwise translate the whole sequences (3 ORFs)
            if 'CDS' in key_values.keys():
                try:
                    cds_info = [int(x) for x in key_values['CDS'].split('-')]
                    ref_seq = ref_seq[cds_info[0]-1:cds_info[1]]
                    this_num_orfs = 1
                    this_num_orfs_complement = 0
                except (ValueError, IndexError, KeyError):
                    print("Could not extra cds position from fasta header for: ", transcript_id, desc) 
            
            ref_orfs = get_orfs(ref_seq, trans_table, this_num_orfs, this_num_orfs_complement)
            
            write_output(seq_id = transcript_id, 
                         description = desc, 
                         seqs= ref_orfs, 
                         prots_fn = prots_fn)
            
    return proteindb_output_file
    
    
def write_output(seq_id, description, seqs, prots_fn):
    "write the orfs to the output file"
    write_i = False
    if len(seqs)>1: #only add _num when multiple ORFs are generated (e.g in 3 ORF)
        write_i = True
    
    for i, orf in enumerate(seqs):
        if write_i:
            prots_fn.write('>{} {}\n{}\n'.format(seq_id+"_"+str(i+1), description, orf))
        else:
            prots_fn.write('>{} {}\n{}\n'.format(seq_id, description, orf))
    
    
if __name__ == '__main__':
    translate_transcript(args.transcripts_fasta,
                  args.proteindb_output_file,
                  args.trans_table, args.num_orfs, args.num_orfs_complement,
                  args.biotype_str,
                  args.excluded_biotypes, args.included_biotypes,
                  args.skip_including_all_CDSs,
                  args.expression_str, args.expression_thresh
                  )
    
    
    