from click.testing import CliRunner

from pypgatk import cli


def vcf_to_proteindb():
    """
    Test the default behaviour of the vep-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli, ['vep-to-proteindb','--config_file', '../config/ensembl_config.yaml', '--vep_annotated_vcf', '../testdata/test.vcf', '--transcript_fasta', '../testdata/test_transcripts.fa', '--gene_annotated_gtf', '../testdata/test.gtf', '--output_proteindb', '../testdata/proteindb.fa'])
    assert result.exit_code == 0

if __name__ == '__main__':
    vcf_to_proteindb()