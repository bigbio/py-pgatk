from click.testing import CliRunner

from pypgatk import cli

def vcf_to_proteindb():
    """
    Test the default behaviour of the vcf-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['vcf-to-proteindb', '--config_file', '../config/ensembl_config.yaml', 
                            '--vep_annotated_vcf', '../testdata/test.vcf', 
                            '--input_fasta', '../testdata/test.fa',
                            '--gene_annotations_gtf', '../testdata/test.gtf', 
                            '--output_proteindb', '../testdata/proteindb_from_ENSEMBL_VCF.fa'])
    assert result.exit_code == 0
    
def vcf_gnomad_to_proteindb():
    """
    Test the default behaviour of the vcf-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['vcf-to-proteindb', '--config_file', '../config/ensembl_config.yaml', 
                            '--vep_annotated_vcf', '../testdata/test_gnomad.vcf', 
                            '--input_fasta', '../testdata/test_gencode.fa',
                            '--gene_annotations_gtf', '../testdata/test_gencode.gtf', 
                            '--output_proteindb', '../testdata/proteindb_from_gnomad_VCF.fa',
                            '--af_field', 'controls_AF', 
                            '--transcript_index', 6, 
                            'biotype_str', 'transcript_type',
                            'annotation_field_name', 'vep'])
    assert result.exit_code == 0

def dnaseq_to_proteindb():
    """
    Test the default behaviour (translate CDSs) of the dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', '../config/ensembl_config.yaml', 
                            '--input_fasta', '../testdata/test.fa', 
                            '--output_proteindb', '../testdata/proteindb_from_CDSs_DNAseq.fa'])
    assert result.exit_code == 0
    
def dnaseq_lncRNAs_to_proteindb():
    """
    Test generation of proteinDB from long noncoding RNAs using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', '../config/ensembl_config.yaml', 
                            '--input_fasta', '../testdata/test.fa', 
                            '--output_proteindb', '../testdata/proteindb_from_lncRNAs_DNAseq.fa',
                            '--include_biotypes', '3prime overlapping ncrna, ambiguous orf, antisense, antisense RNA, lincRNA, ncrna host, processed transcript, sense intronic, sense overlapping', 
                            '--skip_including_all_cds'])
    assert result.exit_code == 0

def dnaseq_sncRNAs_to_proteindb():
    """
    Test generation of proteinDB from short noncoding RNAs using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', '../config/ensembl_config.yaml', 
                            '--input_fasta', '../testdata/test.fa', 
                            '--output_proteindb', '../testdata/proteindb_from_lncRNAs_DNAseq.fa',
                            '--include_biotypes', 'miRNA, miRNA_pseudogene, miscRNA, miscRNA pseudogene, Mt rRNA, Mt tRNA, rRNA, scRNA, snlRNA, snoRNA, snRNA, tRNA, tRNA_pseudogene', 
                            '--skip_including_all_cds'])
    assert result.exit_code == 0

def dnaseq_pseudogenes_to_proteindb():
    """
    Test generation of proteinDB from pseudogenes using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', '../config/ensembl_config.yaml', 
                            '--input_fasta', '../testdata/test.fa', 
                            '--output_proteindb', '../testdata/proteindb_from_pseudogenes_DNAseq.fa',
                            '--include_biotypes', 'disrupted domain, IGC pseudogene, IGJ pseudogene, IG pseudogene, IGV pseudogene, processed pseudogene, transcribed processed pseudogene, transcribed unitary pseudogene, transcribed unprocessed pseudogene, translated processed pseudogene, TRJ pseudogene, unprocessed pseudogene', 
                            '--skip_including_all_cds'])
    assert result.exit_code == 0

def dnaseq_altORFs_to_proteindb():
    """
    Test generation proteinDB from altORFs using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', '../config/ensembl_config.yaml', 
                            '--input_fasta', '../testdata/test.fa', 
                            '--output_proteindb', '../testdata/proteindb_from_altORFs_DNAseq.fa',
                            '--include_biotypes', 'altORFs', '--skip_including_all_cds'])
    assert result.exit_code == 0
    
if __name__ == '__main__':
    
    vcf_to_proteindb()
    vcf_gnomad_to_proteindb
    
    dnaseq_to_proteindb()
    dnaseq_lncRNAs_to_proteindb()
    dnaseq_sncRNAs_to_proteindb()
    dnaseq_pseudogenes_to_proteindb()
    dnaseq_altORFs_to_proteindb()
    