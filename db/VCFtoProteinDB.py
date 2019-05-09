"""
Translate VEP annotated variants from VCF file
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import vcf

def parse_commandline_args():
    
    parser = argparse.ArgumentParser(description='Generate peptides based on DNA variants (VCF)')
    
    #reference annotations and genome files
    parser.add_argument('--genome_fasta', help='path to the genome sequence', required=True)
    parser.add_argument('--gtf_transcripts_fasta', help='path to the transcripts fasta (gtf converted fasta format)', required=True)
    parser.add_argument('--gene_annotations_gtf', help='path to the gene annotations file', required=True)
    
    #input files
    parser.add_argument('--vep_annotated_vcf', help='path to the vep annotated VCF file', required=True)
    
    #translation params
    parser.add_argument('--nuclear_trans_table', default=1, type=int, help='nuclear_trans_table (default 1)')
    parser.add_argument('--mito_trans_table', default=2, type=int, help='mito_trans_table (default 2)')
    
    #output params
    parser.add_argument('--var_prefix', default = "var", help="string to add before the variant peptides")
    parser.add_argument('--report_ref_seq', action='store_const', const=True, help='in addition to var peps, also report all ref peps')
    parser.add_argument('--proteindb_output_file', default = "peps.fa", help="output file name, exits if already exists")
    
    #vcf info fields and filters
    parser.add_argument('--annotation_field_name', default = "CSQ", help="annotation field name found in the INFO column, e.g CSQ or vep")
    parser.add_argument('--AF_field', default = "MAF", help="field name in the VCF INFO column to use for filtering on AF")
    parser.add_argument('--af_threshold', default=0.01, type=float, help='minium AF threshold for considering common variants')
    
    parser.add_argument('--transcript_index', default=3, type=int, help='index of transcript ID in the annotated columns (separated by |)')
    parser.add_argument('--consequence_index', default=1, type=int, help='index of consequence in the annotated columns (separated by |)')
    
    parser.add_argument('--excluded_biotypes', nargs='*', default = [], help="excluded_biotypes")
    parser.add_argument('--excluded_consequences', nargs='*', default = ['downstream_gene_variant', 
                                                                           'upstream_gene_variant', 'intergenic_variant', 
                                                                           'intron_variant', 'synonymous_variant'], 
                                                                           help="excluded_consequences")
    parser.add_argument('--skip_including_all_CDSs', action='store_const', const=True, help="by default any transcript that has a defined CDS will be used, this option disables this features instead it only depends on the biotypes")
    parser.add_argument('--included_biotypes', nargs='*', default = [], 
                        help="included_biotypes")
    parser.add_argument('--included_consequences', nargs='+', default = ['all'], help="included_consequences, default all")
    #GTF params
    parser.add_argument('--biotype_str', default='transcript_biotype', type=str, help='string used to identify gene/transcript biotype in the gtf file.')
    
    return parser.parse_args(sys.argv[1:])
    
    
args = parse_commandline_args()

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


def get_altseq(ref_seq, ref_allele, var_allele, var_pos, strand, features_info, cds_info=[]):
    """
    the given sequence in the fasta file represents all exons of the transcript combined.
    for protein coding genes, the CDS is specified therefore the sequence position has to be 
    calculated based on the CDS's positions
    However, for non-protein coding genes, the whole sequence is used
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
    
    if strand == '-':#get the correct orientation, because exons are oredered based on their position
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
    
    return "", ""
    
    
def parse_gtf(gtf_fn, gtf_db_file):
    """Conver GTF to a FeatureDB"""
    try:
        gffutils.create_db(gtf_fn, gtf_db_file, merge_strategy="create_unique", 
                           keep_order=True, disable_infer_transcripts=True, disable_infer_genes=True,
                           verbose=True,
                           force = False)
    except: #already exists
        pass
    
    db  = gffutils.FeatureDB(gtf_db_file)
    return db
    
    
def get_features(db, feature_id, biotype_str, feature_types = ['exon']):
    """get chr, genomic positions, strand and biotype for feature_id
    also genomic positions for all its elements (exons/cds&start_codon)
    """
    feature = db[feature_id]
    coding_features = []
    features = db.children(feature, featuretype = feature_types, order_by='end')
    for f in features:
        f_type = f.featuretype
        coding_features.append([f.start, f.end, f_type])
    return feature.chrom, feature.strand, coding_features, feature.attributes[biotype_str][0]
    
    
def get_orfs(ref_seq, alt_seq, trans_table, num_orfs=1):
    """translate the coding_ref and the coding_alt into ORFs"""
    
    ref_orfs = []
    alt_orfs = []
    for n in range(0, num_orfs):
        ref_orfs.append(ref_seq[n::].translate(trans_table))
        alt_orfs.append(alt_seq[n::].translate(trans_table))
    
    return ref_orfs, alt_orfs
    
    
def get_protseq_from_vcf(vep_annotated_vcf, 
                      transcripts_fasta,  
                      gtf_db,
                      proteindb_output_file,  
                      nuclear_trans_table, mito_trans_table,
                      AF_field, af_threshold, 
                      annotation_field_name, consequence_index, transcript_index, 
                      biotype_str,
                      excluded_consequences, included_consequences, excluded_biotypes, included_biotypes,
                      skip_including_all_CDSs,
                      var_prefix, report_ref_seq
                      ):
    """
    Generate peps for variants by modifying sequences of affected transcripts (VEP annotated).
    It only considers variants within potential coding regions of the transcript 
    (CDSs & stop codons for protein-coding genes, exons for non-protein coding genes). 
    """
    transcripts_dict = SeqIO.index(transcripts_fasta, "fasta")
    #handle cases where the transript has version in the GTF but not in the VCF
    transcript_id_mapping = {k.split('.')[0]: k for k in transcripts_dict.keys()}    
    
    with open(proteindb_output_file , 'w') as prots_fn:
        vcf_reader = vcf.Reader(open(vep_annotated_vcf, 'r'))
        for record in vcf_reader:
            #only process variants above a given allele frequencey threshold
            #get AF from the INFO field
            try:
                af = float(record.INFO[AF_field])
            except TypeError:
                af = float(record.INFO[AF_field][0])
            except KeyError:
                continue
            
            #check if the AF passed the threshold
            if af < af_threshold:
                continue
            
            trans_table = nuclear_trans_table
            consequences = []
            if str(record.CHROM).lstrip('chr').upper() in ['M', 'MT']:
                trans_table = mito_trans_table
            
            processed_transcript_allele = []
            for transcript_record in record.INFO[annotation_field_name]:
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
                    transcript_id_v = transcript_id_mapping[transcript_id]
                except KeyError:
                    transcript_id_v = transcript_id
                
                try:
                    row = transcripts_dict[transcript_id_v]
                    ref_seq = row.seq#get the seq and desc for the transcript from the fasta of the gtf
                    desc = str(row.description)
                except KeyError:
                    print("Transcript {} not found in fasta of the GTF file {}".format(transcript_id_v, record))
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
                
                chrom, strand, features_info, feature_biotype = get_features(gtf_db, transcript_id_v, biotype_str, feature_types)
                
                #skip transcripts with unwanted consequences
                if (consequence in excluded_consequences or
                    (consequence not in included_consequences and included_consequences != ['all'])
                    ):
                    continue
                
                #only include features that have the specified biotypes or they have CDSs info
                if 'CDS' in feature_types and not skip_including_all_CDSs:
                    pass
                elif (feature_biotype in excluded_biotypes or 
                    (feature_biotype not in included_biotypes and included_biotypes!=['all'])):
                    continue
                
                for alt in record.ALT:#in cases of multiple alternative alleles consider all
                    if transcript_id+str(alt) not in processed_transcript_allele: #because VEP reports affected transcripts per alt allele
                        processed_transcript_allele.append(transcript_id+str(alt))
                        "for non-CDSs, only consider the exon that actually overlaps the variant"
                        if (chrom.lstrip("chr") == str(record.CHROM).lstrip("chr") and
                            check_overlap(record.POS, record.POS+len(alt), features_info)): 
                            coding_ref_seq, coding_alt_seq = get_altseq(ref_seq, Seq(str(record.REF)), Seq(str(alt)), 
                                                 int(record.POS), strand, features_info, cds_info)
                            if coding_alt_seq!="":
                                ref_orfs, alt_orfs = get_orfs(coding_ref_seq, coding_alt_seq, trans_table, num_orfs)
                                
                                record_id = ""
                                if record.ID:
                                    record_id = '_' + str(record.ID)
                                write_output(seq_id = '_'.join([var_prefix + str(record_id), 
                                                       '.'.join([str(record.CHROM), str(record.POS), str(record.REF), str(alt)]), 
                                                       transcript_id_v]), 
                                            desc = feature_biotype+":"+consequence, 
                                            seqs= alt_orfs, 
                                            prots_fn = prots_fn)
                                
                                if report_ref_seq:
                                    write_output(seq_id = transcript_id_v, 
                                                 desc = feature_biotype, 
                                                 seqs= ref_orfs, 
                                                 prots_fn = prots_fn)
    
    return proteindb_output_file


def write_output(seq_id, desc, seqs, prots_fn):
    "write the orfs to the output file"
    for i, alt_orf in enumerate(seqs):
        orf_num = ""
        if i>0:#only add _num when multiple ORFs are generated (e.g in 3 ORF)
            orf_num = "_"+str(i+1)
        
        prots_fn.write('>{} {}\n{}\n'.format(seq_id+orf_num, desc, alt_orf))
    
    
if __name__ == '__main__':
    "generate an sqlite3 db from the GTF file"
    gtf_db  = parse_gtf(args.gene_annotations_gtf, 
                        args.gene_annotations_gtf.replace('.gtf', '.db') )
    "run VCF processor"
    get_protseq_from_vcf(args.vep_annotated_vcf, args.gtf_transcripts_fasta, gtf_db,
                      args.proteindb_output_file,
                      args.nuclear_trans_table, args.mito_trans_table,
                      args.AF_field, args.af_threshold, 
                      args.annotation_field_name, args.consequence_index, args.transcript_index, 
                      args.biotype_str,
                      args.excluded_consequences, args.included_consequences, args.excluded_biotypes, args.included_biotypes,
                      args.skip_including_all_CDSs,
                      args.var_prefix, args.report_ref_seq)
    
    
    