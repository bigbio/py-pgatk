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
                            '--transcript_fasta', '../testdata/test_transcripts.fa',
                            '--gene_annotations_gtf', '../testdata/test.gtf', 
                            '--output_proteindb', '../testdata/proteindb_fromVCF.fa'])
    assert result.exit_code == 0

def dnaseq_to_proteindb():
    """
    Test the default behaviour of the dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', '../config/ensembl_config.yaml', 
                            '--dnaseq_fasta', '../testdata/test_transcripts.fa', 
                            '--output_proteindb', '../testdata/proteindb_fromDNAseq.fa'])
    assert result.exit_code == 0
    
if __name__ == '__main__':
    vcf_to_proteindb()
    dnaseq_to_proteindb()