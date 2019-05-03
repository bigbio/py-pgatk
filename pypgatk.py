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
@click.option('--config-file',
              '-c',
              help='Configuration file for the ensembl data downloader pipeline',
              default='config/ensembl_downloader_config.yaml')
@click.option('--output-directory',
              '-o',
              help='Output directory for the peptide databases',
              default="./database_ensembl/")
@click.option('--folder-prefix-release',
              '-fp', help='Output folder prefix to download the data',
              default='release-')
@click.option('--taxonomy',
              '-t',
              help='Taxonomy List (comma separated) that will be use to download the data from Ensembl',
              default='')
@click.option('--skip-gtf', '-sg', help="Skip the gtf file during the download", is_flag=True)
@click.option('--skip-protein', '-sp', help="Skip the protein fasta file during download", is_flag=True)
@click.option('--skip-cds', '-sc', help='Skip the CDS file download', is_flag=True)
@click.option('--skip-ncrna', '-snr', help='Skip the ncRNA file download', is_flag=True)
@click.pass_context
def ensembl_downloader(ctx, config_file, output_directory, folder_prefix_release, taxonomy, skip_gtf,
                       skip_protein, skip_cds, skip_ncrna):
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
@click.option('--config-file',
              '-c',
              help='Configuration file for the ensembl data downloader pipeline',
              default='config/cbioportal_config.yaml')
@click.option('--output-directory',
              '-o',
              help='Output directory for the peptide databases',
              default="./database_cbioportal/")
@click.option('--list-studies', '-l',
              help='Print the list of all the studies in cBioPortal (https://www.cbioportal.org)', is_flag=True)
@click.option('--download-study', '-d',
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
@click.option('--config-file',
              '-c',
              help='Configuration file for the ensembl data downloader pipeline',
              default='config/cosmic_config.yaml')
@click.option('--output-directory',
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
@click.option('--config-file',
              '-c',
              help='Configuration file for the cosmic data pipelines',
              default='config/cosmic_config.yaml')
@click.option('-in', '--input-mutation', help='Cosmic Mutation data file')
@click.option('-fa', '--input-genes', help='All Cosmic genes')
@click.option('-out', '--output-db', help='Protein database including all the mutations')
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
@click.option('--config-file',
              '-c',
              help='Configuration for ',
              default='config/cbioportal_config.yaml')
@click.option('-in', '--input-mutation', help='Cbioportal mutation file')
@click.option('-fa', '--input-cds', help='CDS genes from ENSEMBL database')
@click.option('-out', '--output-db', help='Protein database including all the mutations')
@click.pass_context
def cbioportal_to_proteindb(ctx, config_file, input_mutation, input_cds, output_db):
    if input_mutation is None or input_cds is None or output_db is None:
        print_help()

    pipeline_arguments = {}
    pipeline_arguments[CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE] = input_mutation
    pipeline_arguments[CancerGenomesService.CONFIG_COMPLETE_GENES_FILE] = input_cds
    pipeline_arguments[CancerGenomesService.CONFIG_OUTPUT_FILE] = output_db

    cosmic_to_proteindb_service = CancerGenomesService(config_file, pipeline_arguments)
    cosmic_to_proteindb_service.cbioportal_to_proteindb()


if __name__ == "__main__":
    cli()
