"""
This is the main tool that give access to all commands and options provided by the py-pgatk

@author ypriverol

"""
import logging
import click
from ensembl.data_downloader import EnsemblDataDownloadService
from toolbox.exceptions import AppConfigException

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


# Cli returns command line requests
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """This is the main tool that give access to all commands and options provided by the pypgatk"""


@cli.command()
@click.option('--config_file',
              '-c',
              help='Configuration file for the ensembl data downloader pipeline',
              default='config/config_ensembl_downloader.yaml')
@click.option('--output_directory',
              '-o',
              help='Output directory for the peptide databases',
              default="./database/")
@click.option('--folder_prefix_release',
              '-fp', help='Output folder prefix to download the data',
              default='release-')
@click.option('--taxonomy',
              '-t',
              help='Taxonomy List (comma separated) that will be use to download the data from Ensembl',
              default='')
@click.pass_context
def ensembl_downloader(ctx, config_file, output_directory, folder_prefix_release, taxonomy):
    """ This tool enables to download from enseml ftp the FASTA and GTF files"""

    if config_file is None:
        msg = "The config file for the pipeline is missing, please provide one "
        logging.error(msg)
        raise AppConfigException(msg)

    pipeline_arguments = {}
    if output_directory is not None:
        pipeline_arguments[EnsemblDataDownloadService._CONFIG_OUTPUT_DIRECTORY] = output_directory
    if folder_prefix_release is not None:
        pipeline_arguments[EnsemblDataDownloadService._CONFIG_KEY_FOLDER_PREFIX_RELEASE] = folder_prefix_release
    if taxonomy is not None:
        pipeline_arguments[EnsemblDataDownloadService._CONFIG_TAXONOMY] = taxonomy

    ensembl_download_service = EnsemblDataDownloadService(config_file, pipeline_arguments)

    logger = ensembl_download_service.get_logger_for("Main Pipeline Ensembl Downloader")
    logger.info("Pipeline STARTING ... ")

    ensembl_download_service.download_database_by_species()

    logger.info("Pipeline Finish !!!")


if __name__ == "__main__":
    cli()
