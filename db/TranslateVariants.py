'''
Translate annotated variants from a given VCF file

-Parse the corresponding GTF
-Convert the GTF file into FASTA
-Store the transcripts' DNA sequences in a dictionary 
-For each variant in the VCF, iteraterate through all listed transcripts
for those that pass the filters, get their DNA sequences
-introduce the alternative allele in the DNA sequence and translate
- Generate the list of peptideas based on the ref seq and the mutated seq
- filter out peps that are found in the ref seq list
'''

import sys, os
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import vcf
import gffutils

    
def check_overlap(var_start, var_end, features_info = [[0, 1, 'type']]):
    "return true when the variant ovelaps any of the features"
    if var_start==-1:
        return True
    #check if the var overlaps any of the features
    for feature_pos in features_info:
        pep_start = feature_pos[0]
        pep_end = feature_pos[1]
        if pep_start>=var_start and pep_start<=var_end:#fully contained or partial overlap from the end
            return True
        elif pep_end>=var_start and pep_end<=var_end:#partial overlap in the begining
            return True
        elif pep_start<=var_start and pep_end>=var_end: #fully covered
            return True
    return False

def parse_commandline_args():
    
    parser = argparse.ArgumentParser(description='Generate peptide databases based on RNA-seq data')
    
    #reference annotations and genome files
    parser.add_argument('--canonical_peps', help='path to the canonical peps sequences', required=True)
    parser.add_argument('--genome_fasta', help='path to the genome sequence', required=True)
    parser.add_argument('--chr_sizes_file', help='path to the chromosome sizes tsv file for the used assembly')
    
    #input files
    parser.add_argument('--vep_annotated_vcf', help='path to the vep annotated VCF file', required=True)
    parser.add_argument('--gene_annotations_gtf', help='path to the gene annotations file')
    
    #exec paths
    parser.add_argument('--stringtie_dir_path', default = "", help="patht to the directory where stringtie and gffread is installed, leave empty if it is available system wide")
    parser.add_argument('--no_exec', action='store_const', const=True, help="Don't run any system commands just print them, for testing purposes")
    
    #vcf processing params
    parser.add_argument('--left_extend_var', default=50, type=int, help='number of bases to subtract before the start of the variant for extracting fasta sequences')
    parser.add_argument('--right_extend_var', default=50, type=int, help='number of bases to add to the variant for extracting fasta sequences')
    
    #translation params
    parser.add_argument('--nuclear_trans_table', help='nuclear_trans_table (default 1)', default=1, type=int)
    parser.add_argument('--mito_trans_table', help='mito_trans_table (default 1)', default=2, type=int)
    
    #digestion params
    parser.add_argument('--min_pep_len', help='minimum tryptic peptide length (default 8)', default=8, type=int)
    parser.add_argument('--max_pep_len', help='maximum tryptic peptide length (default 40)', default=40, type=int)
    
    #output params
    parser.add_argument('--varpep_prefix', default = "var", help="string to add before the variant peptides")
    parser.add_argument('--only_incl_var_peps_as_var', action='store_const', const=True, help='only mark the modified peptides as var and not the ref peps from the region')
    parser.add_argument('--report_ref_peptides', action='store_const', const=True, help='in addition to var peps, also report all ref peps')
    parser.add_argument('--peps_output_file', default = "peps.fa", help="output file name, exits if already exists")
    
    #printing and testing
    parser.add_argument('--verbose', action='store_const', const=True, help='Set to false to avoid pringting commands and file names at each stage of theprocess')
    
    return parser.parse_args(sys.argv[1:])
    
args = parse_commandline_args()

def get_id(record_id):
    return record_id.split('.')[0]

def generate_fasta_from_gtf(gtf_file, fasta_output):
    """convert transcript coordinates to genome fasta sequences"""
    if os.path.isfile(fasta_output):
        return fasta_output
    #add genomic position to the info column so that it is written to the output file
    cmd = "{stringtie_dir_path}gffread -F -w {fasta_out} -g {genome_fasta} {gtf_file}".format(
        stringtie_dir_path=args.stringtie_dir_path, fasta_out=fasta_output, genome_fasta=args.genome_fasta, gtf_file=gtf_file)
    exec(cmd)
    return fasta_output
    
def fasta_to_dict(fasta_file, remove_version_from_id=False):
    genome_dict = {}
    if remove_version_from_id:
        genome_dict = SeqIO.index(fasta_file, "fasta", key_function=get_id)
    else:
        genome_dict = SeqIO.index(fasta_file, "fasta")
    return genome_dict

def get_altseq(ref_seq, ref_allele, var_allele, var_pos, strand, features_info, cds_info=[]):
    """
    the given sequence in the fasta file represents all exons of the transcript combined.
    for protein coding genes, the CDS is specified therefore the sequence position has to be 
    calculated based on the CDS's positions
    However, for non-protein coding genes, the whole sequence should be used
    1. find which exon or cds is affected by the variant.
    deletions??
    """
    alt_seq = ""
    if len(cds_info)==2:
        start_coding_index = cds_info[0]-1 #it should be index not pos
        stop_coding_index = cds_info[1]#get end position of the  last cds
    else:
        start_coding_index = 0
        total_len = 0
        for x in features_info:
            total_len+= x[1]-x[0]+1 
        stop_coding_index = total_len#the features are sorted by end therefroe the end pos of the last item is the last coding nc
    
    if strand == '-':#ge the correct orientation, because exons are oredered based on their position
        ref_seq = ref_seq[::-1] #in order to calculate from the first base of the first feature (sorted by genomic coordinates)
        ref_allele = ref_allele.complement() #the reverse will be done on return
        var_allele = var_allele.complement() #the reverse will be done on return
    
    ref_seq = ref_seq[start_coding_index:stop_coding_index]#just keep the coding regions (mostly effective in case of protein-coding genes) 
    nc_index = 0
    feature_len = 0
    if (len(ref_allele) == len(var_allele) or ref_allele[0]==var_allele[0]):
        for feature in features_info:#for every exon, cds or stop codon
            if var_pos in range(feature[0], feature[1]+1):#get index of the var relative to the position of the overlapping feature in the coding region 
                var_index_in_cds =  nc_index + (var_pos - feature[0])
                #modify the coding reference sequence accoding to the var_allele
                c = len(ref_allele)
                alt_seq = ref_seq[0:var_index_in_cds] + var_allele + ref_seq[var_index_in_cds+c::]#variant and ref strand??
                if strand == '-':
                    return ref_seq[::-1], alt_seq[::-1]
                else:
                    return ref_seq, alt_seq
                
            feature_len = (feature[1]-feature[0]+1)
            nc_index+=feature_len
            #feature_seq = coding_ref_seq[nc_index:nc_index+feature_len]
            
    return ref_seq, alt_seq
    
def parse_gtf(gtf_fn):
    gtf_db_file = gtf_fn.replace('.gtf', '.db') 
    if not os.path.isfile(gtf_db_file):
        gffutils.create_db(gtf_fn, gtf_db, merge_strategy="create_unique", 
                           keep_order=True, disable_infer_transcripts=True, disable_infer_genes=True,
                           verbose=True,)
    db  = gffutils.FeatureDB(gtf_db_file)
    return db
    #trans = db.children(feature_id, featuretype=feature_type)
    
def get_features(db, feature_id, feature_types = ['exon']):
    """returns chromsome and strand for the transcript and 
    all exons/cds&start_codon in the third element of a list
    """
    feature = db[feature_id]
    coding_features = []
    features = db.children(feature, featuretype = feature_types, order_by='end')
    for f in features:
        f_type = f.featuretype
        coding_features.append([f.start, f.end, f_type])
    return feature.chrom, feature.strand, coding_features, feature.attributes['gene_biotype'][0]
    
    
def get_orfs(ref_seq, alt_seq, trans_table, num_orfs=1, num_orfs_complement=0):
    """
    translate and digest the coding_ref and the coding_alt
    remove parts of the coding_alt that are in the coding_ref and return the alt peptides
    """
    ref_orfs = []
    alt_orfs = []
    for n in range(0, num_orfs):
        ref_orfs.append(ref_seq[n::].translate(trans_table))
        alt_orfs.append(alt_seq[n::].translate(trans_table))
    
    rev_ref_seq = ref_seq.reverse_complement()
    rev_alt_seq = alt_seq.reverse_complement()
    for n in range(0, num_orfs_complement):
        ref_orfs.append(rev_ref_seq[n::].translate(trans_table))
        alt_orfs.append(rev_alt_seq[n::].translate(trans_table))
    
    return ref_orfs, alt_orfs
    

def digest_seq(seq, digestion_enzyme = 'trypsin', cut_sites = ['K','R'], miss_cleavage = 0):
    peptide = ''
    peptides_pos = {}
    peptides = []
    
    for c, aa in enumerate(seq):
        peptide+=aa
        next_aa = ''
        try:
            next_aa = seq[c+1]
        except IndexError:
            pass
        if digestion_enzyme=='trypsin': 
            if (aa in cut_sites and next_aa!='P'):
                if len(peptide)>0:
                    if len(peptides) == 0:
                        peptides_pos[peptide] = [0, c*3]
                    else:
                        peptides_pos[peptide] = [(c*3-(len(peptide)*3)-1), c*3]
                    peptides.append(peptide)
                peptide = ''
                continue
        else:
            """No digestion enzyme is used. Peptides between stop codons are generated from the next statement. 
                implement other digestion enzymes here"""
            pass
    if len(peptide)>0:
        peptides.append(peptide)
        if len(peptides) == 0:
            peptides_pos[peptide] = [0, c*3]
        else:
            peptides_pos[peptide] = [(c*3-(len(peptide)*3)-1), c*3]
        
        
    peptides_with_miss_cleavage = []
    for i in range(1, miss_cleavage+1):
        for j,pep in enumerate(peptides):
            if j+i<len(peptides):
                peptide = ''.join([x for x in (peptides[j:j+i+1])])
                peptides_with_miss_cleavage.append(peptide)
                peptides_pos[peptide] = peptides_pos[peptides[j+i+1]]#to be checked
    
    peptides.extend(peptides_with_miss_cleavage)
    
    return peptides, peptides_pos


def get_peps(ref_orfs, alt_orfs, alt,
             var_prefix,
             min_pep_len, max_pep_len,
             only_incl_var_peps_as_var,
             record_id, desc,
             left_extend=0, right_extend=0):
    """
    digest each orf
    """
    all_ref_peptides = []#peps from all ref orfs
    ref_peps_to_return = []
    var_peps_to_return = []
    #get tryptic peptides from the ref orfs
    for i, prot_seq in enumerate(ref_orfs):
        ref_peps = []
        for ref_seq in prot_seq.split('*'):
            peptides, peptides_pos = digest_seq(ref_seq, miss_cleavage = 0)  
            ref_peps.extend(peptides)
            
        for n, ref_pep in enumerate(list(set(ref_peps))):
            if len(ref_pep) in range(min_pep_len, max_pep_len+1) and set(ref_pep)!={'X'}:
                ref_peps_to_return.append('>{}_{}_{} {}\n{}\n'.format(record_id, i+1, n+1, desc, ref_pep))
                all_ref_peptides.append(ref_pep)
    
    #get tryptic peptides from the alt orfs and filter against peptides in the ref set
    for i, prot_seq in enumerate(alt_orfs):
        alt_peps = []
        for alt_seq in prot_seq.split('*'):
            peptides, peptides_pos = digest_seq(alt_seq, miss_cleavage = 0)
            if (left_extend>0 or right_extend>0):
                print(alt_seq)
                for pep in peptides:
                    print(pep)
                    #check the coords
                    if check_overlap(var_start = left_extend-1, 
                                     var_end = left_extend+len(alt), 
                                     features_info = [[
                                         peptides_pos[pep][0], 
                                         peptides_pos[pep][1], 'none']]):
                        print("yesssssss*****", pep)
                        """if the peptide was originated from the variant sequence then include 
                        otherwise ignore the extended part that does not overlap the variant"""
                        alt_peps.append(pep)
            else:
                alt_peps.extend(peptides)
        uniq_alt_peps = set(alt_peps)
        if only_incl_var_peps_as_var:
            uniq_alt_peps = uniq_alt_peps - set(all_ref_peptides)#don't include peptides that are found in the ref orfs
        
        for n, alt_pep in enumerate(uniq_alt_peps):
            if len(alt_pep) in range(min_pep_len, max_pep_len+1) and set(alt_pep)!={'X'}:
                var_peps_to_return.append('>{}{}_{}_{} {}\n{}\n'.format(var_prefix, record_id, i+1, n+1, desc, alt_pep))
                
    return ref_peps_to_return, var_peps_to_return    
    
def write_out_peps(peps_fn, all_ref_peps, all_var_peps):
    
    if args.report_ref_peptides:
        for pep in all_ref_peps:
            peps_fn.write(pep)
    for pep in all_var_peps:
        peps_fn.write(pep)
    
    return [], []
    

def get_peps_from_vcf(vep_annotated_vcf, #VCF input file 
                      transcripts_dict, #dictionary containing transcriptID: seq, desc 
                      gtf_db,   #gffutils GTF sqlite db file
                      transcript_index = 3, #index of transcript ID in the CSQ field of the INFO column
                      consequence_index = 1, #index of consequence in the CSQ field of the INFO column
                      excluded_biotypes = [],
                      excluded_consequences = ['downstream_gene_variant', 'upstream_gene_variant', 
                                               'intergenic_variant', 'intron_variant', 'synonymous_variant'],
                      included_biotypes = ['protein_coding'],
                      included_consequences = ['all'],
                      maf_threshold = 0.01
                      ):
    """
    Generate peps for variants by modifying sequences of affected transcripts (VEP annotated).
    It only considers variants within potential coding regions of the transcript 
    (CDSs & stop codons for protein-coding genes, exons for non-protein coding genes). 
    """
    with open(args.peps_output_file , 'w') as peps_fn:
        all_ref_peps = []
        all_var_peps = []
        vcf_reader = vcf.Reader(open(vep_annotated_vcf, 'r'))
        for r, record in enumerate(vcf_reader):
            #only process variants above a given allele frequencey threshold
            try:
                if record.MAF<maf_threshold:
                    continue
            except KeyError:
                continue
            
            trans_table = args.nuclear_trans_table
            consequences = []
            biotypes = []
            if str(record.CHROM).lstrip('chr').upper() in ['M', 'MT']:
                trans_table = args.mito_trans_table
            
            processed_transcript_allele = []
            for transcript_record in record.INFO['CSQ']:
                transcript_info = transcript_record.split('|')
                
                try:
                    consequence = transcript_info[consequence_index]
                except IndexError:
                    print("Give a valid index for the consequence in the INFO field for: ", transcript_record)
                    continue
                consequences.append(consequence)
                
                try:
                    transcript_id = transcript_info[transcript_index]
                except IndexError:
                    print("Give a valid index for the Transcript ID in the INFO field for: ", transcript_record)
                    continue
                if transcript_id=="":
                    continue
                
                try:
                    row = transcripts_dict[transcript_id]
                    ref_seq = row.seq#get the seq and desc for the transcript from the fasta of the gtf
                    desc = str(row.description)
                except KeyError:
                    print("Transcript {} not found in fasta of the GTF file {}".format(transcript_id, record))
                    continue
                
                feature_types = ['exon']
                #check if cds info exists in the fasta header otherwise translate all exons
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
                
                chrom, strand, features_info, feature_biotype = get_features(gtf_db, transcript_id, feature_types)
                biotypes.append(feature_biotype)
                
                #skip transcripts with unwanted consequences
                if (consequence in excluded_consequences or
                    (consequence not in included_consequences and included_consequences != ['all'])
                    ):
                    continue
                
                #only include features that have the specified biotypes
                if (feature_biotype in excluded_biotypes or 
                    (feature_biotype not in included_biotypes and included_biotypes!=['all'])):
                    continue
                
                for alt in record.ALT:#in cases of multiple alternative alleles consider all
                    if transcript_id+str(alt) not in processed_transcript_allele: #because VEP reports affected transcripts per alt allele
                        processed_transcript_allele.append(transcript_id+str(alt))
                        "for non-CDSs, only consider the exon that actually overlaps the variant"
                        if (chrom.lstrip("chr") ==str(record.CHROM).lstrip("chr") and
                            check_overlap(record.POS, record.POS+len(alt), features_info)): 
                            coding_ref_seq, coding_alt_seq = get_altseq(ref_seq, Seq(str(record.REF)), Seq(str(alt)), 
                                                 int(record.POS), strand, features_info, cds_info)
                            if coding_alt_seq!="":
                                ref_orfs, alt_orfs = get_orfs(coding_ref_seq, coding_alt_seq, trans_table, num_orfs, num_orfs_complement=0)
                                
                                ref_peps, var_peps = get_peps(ref_orfs, alt_orfs, str(alt), args.varpep_prefix, 
                                                args.min_pep_len, args.max_pep_len, args.only_incl_var_peps_as_var,
                                                transcript_id, desc = str(record.ID) + ":" + str(record.POS)+":"+str(alt)+":"+feature_biotype+":"+consequence) 
                                
                                if len(var_peps)>0:
                                    all_var_peps.extend(var_peps)
                                    all_ref_peps.extend(ref_peps)
            
            #write to output periodically
            if len(all_ref_peps)>1000000 or len(all_var_peps)>1000000:
                print('Processed: ', r, ' variant records.')
                all_ref_peps, all_var_peps = write_out_peps(peps_fn, all_ref_peps, all_var_peps)
        #finally write remaining peptides, if any
        all_ref_peps, all_var_peps = write_out_peps(peps_fn, all_ref_peps, all_var_peps)
    
    return args.peps_output_file


if __name__ == '__main__':
    
    gtf_db  = parse_gtf(args.gene_annotations_gtf)
    transcripts_fasta = generate_fasta_from_gtf(args.gene_annotations_gtf, args.gene_annotations_gtf.replace('.gtf', '.fasta'))
    gtf_transcripts_dict = fasta_to_dict(transcripts_fasta)
    get_peps_from_vcf(args.vep_annotated_vcf, gtf_transcripts_dict, gtf_db)
    
    
    
    