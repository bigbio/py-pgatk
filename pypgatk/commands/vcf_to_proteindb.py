import click

from pypgatk.commands.utils import print_help
from pypgatk.ensembl.ensembl import EnsemblDataService


@click.command('vcf-to-proteindb', short_help="Generate peptides based on DNA variants from ENSEMBL VEP VCF files")
@click.option('--config_file', '-c', help='Configuration to perform conversion between ENSEMBL Files',
              default='config/ensembl_config.yaml')
@click.option('--input_fasta', help='Path to the transcript sequence')
@click.option('--vep_annotated_vcf', help='Path to the vep annotated VCF file')
@click.option('--gene_annotations_gtf', help='Path to the gene annotations file')
@click.option('--translation_table', default=1, type=int, help="Translation table (Default 1) ")
@click.option('--mito_translation_table', default=2, type=int, help='Mito_trans_table (default 2)')
@click.option('--var_prefix', default="var", help="String to add before the variant peptides")
@click.option('--report_ref_seq', help='In addition to var peps, also report all ref peps', is_flag=True)
@click.option('--output_proteindb', default="peptide-database.fa", help="Output file name, exits if already exists")
@click.option('--annotation_field_name', default="CSQ",
              help="Annotation field name found in the INFO column, e.g CSQ or vep")
@click.option('--af_field', default="MAF", help="field name in the VCF INFO column to use for filtering on AF")
@click.option('--af_threshold', default=0.01, help='Minium AF threshold for considering common variants')
@click.option('--transcript_index', default=3, type=int,
              help='Index of transcript ID in the annotated columns (separated by |)')
@click.option('--consequence_index', default=1, type=int,
              help='Index of consequence in the annotated columns (separated by |)')
@click.option('--exclude_biotypes', default='', help="Excluded Biotypes")
@click.option('--exclude_consequences',
              default='downstream_gene_variant, upstream_gene_variant, intergenic_variant, intron_variant, synonymous_variant',
              help="Excluded Consequences", show_default=True)
@click.option('--skip_including_all_cds',
              help="by default any transcript that has a defined CDS will be used, this option disables this features instead it only depends on the biotypes",
              is_flag=True)
@click.option('--include_biotypes', default='', help="Include Biotypes")
@click.option('--include_consequences', default='all', help="included_consequences, default all")
@click.option('--biotype_str', default='transcript_biotype',
              help='string used to identify gene/transcript biotype in the gtf file.')
@click.option('--ignore_filters',
              help="enabling this option causes or variants to be parsed. By default only variants that have not failed any filters will be processed (FILTER column is PASS, None, .) or if the filters are subset of the accepted filters. (default is False)",
              is_flag=True)
@click.option('--accepted_filters', default='', help="Accepted filters for variant parsing")
@click.pass_context
def vcf_to_proteindb(ctx, config_file, input_fasta, vep_annotated_vcf, gene_annotations_gtf, translation_table,
                     mito_translation_table,
                     var_prefix, report_ref_seq, output_proteindb, annotation_field_name,
                     af_field, af_threshold, transcript_index, consequence_index, exclude_biotypes,
                     exclude_consequences, skip_including_all_cds, include_biotypes, include_consequences, biotype_str,
                     ignore_filters, accepted_filters):
    if input_fasta is None or vep_annotated_vcf is None or gene_annotations_gtf is None:
        print_help()

    pipeline_arguments = {}
    pipeline_arguments[EnsemblDataService.MITO_TRANSLATION_TABLE] = mito_translation_table
    pipeline_arguments[EnsemblDataService.TRANSLATION_TABLE] = translation_table
    pipeline_arguments[EnsemblDataService.HEADER_VAR_PREFIX] = var_prefix
    pipeline_arguments[EnsemblDataService.REPORT_REFERENCE_SEQ] = report_ref_seq
    pipeline_arguments[EnsemblDataService.PROTEIN_DB_OUTPUT] = output_proteindb
    pipeline_arguments[EnsemblDataService.ANNOTATION_FIELD_NAME] = annotation_field_name
    pipeline_arguments[EnsemblDataService.AF_FIELD] = af_field
    pipeline_arguments[EnsemblDataService.AF_THRESHOLD] = af_threshold
    pipeline_arguments[EnsemblDataService.TRANSCRIPT_INDEX] = transcript_index
    pipeline_arguments[EnsemblDataService.CONSEQUENCE_INDEX] = consequence_index
    pipeline_arguments[EnsemblDataService.EXCLUDE_BIOTYPES] = exclude_biotypes
    pipeline_arguments[EnsemblDataService.EXCLUDE_CONSEQUENCES] = exclude_consequences
    pipeline_arguments[EnsemblDataService.SKIP_INCLUDING_ALL_CDS] = skip_including_all_cds
    pipeline_arguments[EnsemblDataService.INCLUDE_BIOTYPES] = include_biotypes
    pipeline_arguments[EnsemblDataService.INCLUDE_CONSEQUENCES] = include_consequences
    pipeline_arguments[EnsemblDataService.BIOTYPE_STR] = biotype_str
    pipeline_arguments[EnsemblDataService.IGNORE_FILTERS] = ignore_filters
    pipeline_arguments[EnsemblDataService.ACCEPTED_FILTERS] = accepted_filters

    ensembl_data_service = EnsemblDataService(config_file, pipeline_arguments)
    ensembl_data_service.vcf_to_proteindb(vep_annotated_vcf, input_fasta, gene_annotations_gtf)
