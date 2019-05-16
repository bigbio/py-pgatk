"""
This is the main tool that give access to all commands and options provided by the py-pgatk

@author ypriverol

"""
import logging
import click

from cgenomes.cbioportal_downloader import CbioPortalDownloadService
from cgenomes.cgenomes_proteindb import CancerGenomesService
from cgenomes.cosmic_downloader import CosmicDownloadService
from ensembl.data_downloader import EnsemblDataDownloadService
from ensembl.ensembl import EnsemblDataService
from toolbox.exceptions import AppConfigException

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def print_help():
    """
    Print the help of the tool
    :return:
    """
    ctx = click.get_current_context()
    click.echo(ctx.get_help())
    ctx.exit()


# Cli returns command line requests
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """This is the main tool that give access to all commands and options provided by the pypgatk"""


@cli.command('ensembl-downloader', short_help='Command to download the ensembl information')
@click.option('--config_file',
              '-c',
              help='Configuration file for the ensembl data downloader pipeline',
              default='config/ensembl_downloader_config.yaml')
@click.option('--output_directory',
              '-o',
              help='Output directory for the peptide databases',
              default="./database_ensembl/")
@click.option('--folder_prefix_release',
              '-fp', help='Output folder prefix to download the data',
              default='release-')
@click.option('--taxonomy',
              '-t',
              help='Taxonomy List (comma separated) that will be use to download the data from Ensembl',
              default='')
@click.option('--skip_gtf', '-sg', help="Skip the gtf file during the download", is_flag=True)
@click.option('--skip_protein', '-sp', help="Skip the protein fasta file during download", is_flag=True)
@click.option('--skip_cds', '-sc', help='Skip the CDS file download', is_flag=True)
@click.option('--skip_ncrna', '-snr', help='Skip the ncRNA file download', is_flag=True)
@click.option('--skip_vcf', '-vcf', help='Skip the VCF variant file', is_flag=True)
@click.pass_context
def ensembl_downloader(ctx, config_file, output_directory, folder_prefix_release, taxonomy, skip_gtf,
                       skip_protein, skip_cds, skip_ncrna, skip_vcf):
    """ This tool enables to download from enseml ftp the FASTA and GTF files"""

    if config_file is None:
        msg = "The config file for the pipeline is missing, please provide one "
        logging.error(msg)
        raise AppConfigException(msg)

    # Parse pipelines parameters.
    pipeline_arguments = {}
    if output_directory is not None:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_OUTPUT_DIRECTORY] = output_directory
    if folder_prefix_release is not None:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_FOLDER_PREFIX_RELEASE] = folder_prefix_release
    if taxonomy is not None:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_TAXONOMY] = taxonomy
    if skip_protein is not None and skip_protein:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_PROTEIN] = True
    else:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_PROTEIN] = False
    if skip_gtf is not None and skip_gtf:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_GTF] = True
    else:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_GTF] = False
    if skip_cds is not None and skip_cds:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_CDS] = True
    else:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_CDS] = False

    if skip_ncrna is not None and skip_ncrna:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_NCRNA] = True
    else:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_NCRNA] = False

    ensembl_download_service = EnsemblDataDownloadService(config_file, pipeline_arguments)

    logger = ensembl_download_service.get_logger_for("Main Pipeline Ensembl Downloader")
    logger.info("Pipeline STARTING ... ")

    ensembl_download_service.download_database_by_species()

    logger.info("Pipeline Finish !!!")


@cli.command('cbioportal-downloader', short_help=' Command to download the the cbioportal studies')
@click.option('--config_file', '-c', help='Configuration file for the ensembl data downloader pipeline',
              default='config/cbioportal_config.yaml')
@click.option('--output_directory', '-o', help='Output directory for the peptide databases',
              default="./database_cbioportal/")
@click.option('--list_studies', '-l',
              help='Print the list of all the studies in cBioPortal (https://www.cbioportal.org)', is_flag=True)
@click.option('--download_study', '-d',
              help="Download an specific Study from cBioPortal -- (all to download all studies)")
@click.pass_context
def cbioportal_downloader(ctx, config_file, output_directory, list_studies, download_study):
    if config_file is None:
        msg = "The config file for the pipeline is missing, please provide one "
        logging.error(msg)
        raise AppConfigException(msg)

    pipeline_arguments = {}
    if output_directory is not None:
        pipeline_arguments[CbioPortalDownloadService.CONFIG_OUTPUT_DIRECTORY] = output_directory
    if list_studies is not None:
        pipeline_arguments[CbioPortalDownloadService.CONFIG_LIST_STUDIES] = list_studies

    cbioportal_downloader_service = CbioPortalDownloadService(config_file, pipeline_arguments)

    if list_studies is not None:
        list_studies = cbioportal_downloader_service.get_cancer_studies()
        print(list_studies)

    if download_study is not None:
        cbioportal_downloader_service.download_study(download_study)


@cli.command('cosmic-downloader', short_help='Command to download the cosmic mutation database')
@click.option('--config_file',
              '-c',
              help='Configuration file for the ensembl data downloader pipeline',
              default='config/cosmic_config.yaml')
@click.option('--output_directory',
              '-o',
              help='Output directory for the peptide databases',
              default="./database_cosmic/")
@click.option('--username', '-u',
              help="Username for cosmic database -- please if you don't have one register here (https://cancer.sanger.ac.uk/cosmic/register)")
@click.option('--password', '-p',
              help="Password for cosmic database -- please if you don't have one register here (https://cancer.sanger.ac.uk/cosmic/register)")
@click.pass_context
def cosmic_downloader(ctx, config_file, output_directory, username, password):
    if config_file is None:
        msg = "The config file for the pipeline is missing, please provide one "
        logging.error(msg)
        raise AppConfigException(msg)

    pipeline_arguments = {}
    if output_directory is not None:
        pipeline_arguments[CosmicDownloadService.CONFIG_OUTPUT_DIRECTORY] = output_directory
    if username is not None:
        pipeline_arguments[CosmicDownloadService.CONFIG_COSMIC_FTP_USER] = username
    else:
        print_help()

    if password is not None:
        pipeline_arguments[CosmicDownloadService.CONFIG_COSMIC_FTP_PASSWORD] = password
    else:
        print_help()

    cosmic_downloader_service = CosmicDownloadService(config_file, pipeline_arguments)

    cosmic_downloader_service.download_mutation_file()


@cli.command('cosmic-to-proteindb', short_help='Command to translate Cosmic mutation data into proteindb')
@click.option('--config_file',
              '-c',
              help='Configuration file for the cosmic data pipelines',
              default='config/cosmic_config.yaml')
@click.option('-in', '--input_mutation', help='Cosmic Mutation data file')
@click.option('-fa', '--input_genes', help='All Cosmic genes')
@click.option('-out', '--output_db', help='Protein database including all the mutations')
@click.pass_context
def cosmic_to_proteindb(ctx, config_file, input_mutation, input_genes, output_db):
    if input_mutation is None or input_genes is None or output_db is None:
        print_help()

    pipeline_arguments = {CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE: input_mutation,
                          CancerGenomesService.CONFIG_COMPLETE_GENES_FILE: input_genes,
                          CancerGenomesService.CONFIG_OUTPUT_FILE: output_db}

    cosmic_to_proteindb_service = CancerGenomesService(config_file, pipeline_arguments)
    cosmic_to_proteindb_service.cosmic_to_proteindb()


@cli.command('cbioportal-to-proteindb', short_help='Command to translate cbioportal mutation data into proteindb')
@click.option('--config_file',
              '-c',
              help='Configuration for cbioportal to proteindb tool',
              default='config/cbioportal_config.yaml')
@click.option('-in', '--input_mutation', help='Cbioportal mutation file')
@click.option('-fa', '--input_cds', help='CDS genes from ENSEMBL database')
@click.option('-out', '--output_db', help='Protein database including all the mutations')
@click.pass_context
def cbioportal_to_proteindb(ctx, config_file, input_mutation, input_cds, output_db):
    if input_mutation is None or input_cds is None or output_db is None:
        print_help()

    pipeline_arguments = {CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE: input_mutation,
                          CancerGenomesService.CONFIG_COMPLETE_GENES_FILE: input_cds,
                          CancerGenomesService.CONFIG_OUTPUT_FILE: output_db}

    cosmic_to_proteindb_service = CancerGenomesService(config_file, pipeline_arguments)
    cosmic_to_proteindb_service.cbioportal_to_proteindb()


@cli.command('threeframe-translation', short_help="Command to perform 3'frame translation")
@click.option('--config_file',
              '-c',
              help='Configuration to perform conversion between ENSEMBL Files',
              default='config/ensembl_config.yaml')
@click.option('-in', '--input', help='input file to perform the translation')
@click.option('-t', '--translation_table', help='Translation table default value 1', default='1')
@click.option('-out', '--output', help='Output File', default="peptide-database.fa")
@click.pass_context
def threeframe_translation(ctx, config_file, input, translation_table, output):
    if input is None:
        print_help()

    pipeline_arguments = {EnsemblDataService.TRANSLATION_TABLE: translation_table,
                          EnsemblDataService.PROTEIN_DB_OUTPUT: output}

    ensembl_data_service = EnsemblDataService(config_file, pipeline_arguments)
    ensembl_data_service.three_frame_translation(input)


@cli.command('vcf-to-proteindb', short_help="Generate peptides based on DNA variants from ENSEMBL VEP VCF files")
@click.option('--config_file', '-c', help='Configuration to perform conversion between ENSEMBL Files',
              default='config/ensembl_config.yaml')
@click.option('--transcripts_fasta', help='Path to the transcript sequence')
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
@click.pass_context
def vcf_to_proteindb(ctx, config_file, transcripts_fasta, vep_annotated_vcf, gene_annotations_gtf, translation_table,
                     mito_translation_table,
                     var_prefix, report_ref_seq, output_proteindb, annotation_field_name,
                     af_field, af_threshold, transcript_index, consequence_index, exclude_biotypes,
                     exclude_consequences, skip_including_all_cds, include_biotypes, include_consequences, biotype_str):
    if transcripts_fasta is None or vep_annotated_vcf is None or gene_annotations_gtf is None:
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

    ensembl_data_service = EnsemblDataService(config_file, pipeline_arguments)
    ensembl_data_service.vep_to_proteindb(vep_annotated_vcf, transcripts_fasta, gene_annotations_gtf)


@cli.command("dnaseq-to-proteindb", short_help="Generate peptides based on DNA sequences")
@click.option('--config_file', '-c', help='Configuration to perform conversion between ENSEMBL Files',
              default='config/ensembl_config.yaml')
@click.option('--dnaseq_fasta', help='Path to sequences fasta')
@click.option('--translation_table', default=1, type=int, help='Translation Table (default 1)')
@click.option('--num_orfs', default=3, type=int, help='Number of ORFs (default 0)')
@click.option('--num_orfs_complement', default=0, type=int,
              help='Number of ORFs from the reverse side (default 0)')
@click.option('--output_proteindb', default="peptide-database.fa", help="Output file name, exits if already exists")
@click.option('--skip_including_all_cds',
              help="By default any transcript that has a defined CDS will be translated, this option disables this features instead it only depends on the biotypes",
              is_flag=True)
@click.option('--include_biotypes', default='', help="Include Biotypes")
@click.option('--exclude_biotypes', default='', help="Exclude Biotypes")
@click.option('--biotype_str', default='transcript_biotype', type=str,
              help='String used to identify gene/transcript biotype in the gtf file.')
@click.option('--expression_str', default="", type=str,
              help='String to be used for extracting expression value (TPM, FPKM, etc).')
@click.option('--expression_thresh', default=5.0, type=float,
              help='Threshold used to filter transcripts based on their expression values')
@click.pass_context
def dnaseq_to_proteindb(ctx, config_file, dnaseq_fasta, translation_table, num_orfs, num_orfs_complement,
                     output_proteindb,
                     skip_including_all_cds, include_biotypes, exclude_biotypes, biotype_str, expression_str,
                     expression_thresh):
    if dnaseq_fasta is None:
        print_help()

    pipeline_arguments = {}
    pipeline_arguments[EnsemblDataService.TRANSLATION_TABLE] = translation_table
    pipeline_arguments[EnsemblDataService.PROTEIN_DB_OUTPUT] = output_proteindb
    pipeline_arguments[EnsemblDataService.EXCLUDE_BIOTYPES] = exclude_biotypes
    pipeline_arguments[EnsemblDataService.SKIP_INCLUDING_ALL_CDS] = skip_including_all_cds
    pipeline_arguments[EnsemblDataService.INCLUDE_BIOTYPES] = include_biotypes
    pipeline_arguments[EnsemblDataService.BIOTYPE_STR] = biotype_str
    pipeline_arguments[EnsemblDataService.NUM_ORFS] = num_orfs
    pipeline_arguments[EnsemblDataService.NUM_ORFS_COMPLEMENT] = num_orfs_complement
    pipeline_arguments[EnsemblDataService.EXPRESSION_STR] = expression_str
    pipeline_arguments[EnsemblDataService.EXPRESSION_THRESH] = expression_thresh

    ensembl_data_service = EnsemblDataService(config_file, pipeline_arguments)
    ensembl_data_service.dnaseq_to_proteindb(dnaseq_fasta)


if __name__ == "__main__":
    cli()
