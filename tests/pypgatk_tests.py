from click.testing import CliRunner

from pypgatk import cli


# Tests for vcf-to-proteindb
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
                            '--transcript_fasta', '../testdata/test_transcripts_gencode.fa',
                            '--gene_annotations_gtf', '../testdata/test_gencode.gtf',
                            '--output_proteindb', '../testdata/proteindb_from_gnomad_VCF.fa',
                            '--af_field', 'AF',
                            '--transcript_index', 6,
                            'biotype_str', 'transcript_type',
                            'annotation_field_name', 'vep'])
    assert result.exit_code == 0


# Tests for dnaseq-to-proteindb
def dnaseq_to_proteindb():
    """
    Test the default behaviour (translate CDSs) of the dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', '../config/ensembl_config.yaml',
                            '--dnaseq_fasta', '../testdata/test_transcripts.fa',
                            '--output_proteindb', '../testdata/proteindb_from_CDSs_DNAseq.fa'])
    assert result.exit_code == 0


def dnaseq_altorfs_to_proteindb():
    """
    Test generation proteinDB from altORFs using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', '../config/ensembl_config.yaml',
                            '--dnaseq_fasta', '../testdata/test_transcripts.fa',
                            '--output_proteindb', '../testdata/proteindb_from_altORFs_DNAseq.fa',
                            '--include_biotypes', 'altORFs', '--skip_including_all_cds'])
    assert result.exit_code == 0


def dnaseq_lncRNAs_to_proteindb():
    """
    Test generation of proteinDB from lncRNAs using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', '../config/ensembl_config.yaml',
                            '--dnaseq_fasta', '../testdata/test_transcripts.fa',
                            '--output_proteindb', '../testdata/proteindb_from_lncRNAs_DNAseq.fa',
                            '--include_biotypes', 'lncRNA', '--skip_including_all_cds'])
    assert result.exit_code == 0


def dnaseq_pseudogenes_to_proteindb():
    """
    Test generation of proteinDB from pseudogenes using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', '../config/ensembl_config.yaml',
                            '--dnaseq_fasta', '../testdata/test_transcripts.fa',
                            '--output_proteindb', '../testdata/proteindb_from_pseudogenes_DNAseq.fa',
                            '--include_biotypes', 'pseudogenes', '--skip_including_all_cds'])
    assert result.exit_code == 0


if __name__ == '__main__':
    vcf_to_proteindb()
    dnaseq_to_proteindb()
