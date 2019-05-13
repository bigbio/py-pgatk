from click.testing import CliRunner

from pypgatk import cli


def vcf_to_proteindb():

    runner = CliRunner()
    result = runner.invoke(cli, ['--vep_annotated_vcf', '../testdata/test.vcf', '--genome_fasta', '../testdata/test_transcripts.fa', '--gene_annotated_gtf', '../testdata/test.gft'])
    assert result.exit_code == 0
    assert result.output == 'Opt: An Option  Arg: An Arg\n'


if __name__ == '__main__':
    vcf_to_proteindb()