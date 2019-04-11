"""
This pipeline collects data from Ensembl for proteogenomics studies. The pipeline has multiple configurable
options

@author ypriverol

"""
import logging
import click
from ensembl.data_downloader import EnsemblDataDownloadService
from toolbox.exceptions import AppConfigException

__configuration_file = None
__pipeline_arguments = None
__pipeline_director = None


def print_help(ctx, param, value):
    if value is False:
        return
    click.echo(ctx.get_help())
    ctx.exit()


def get_pipeline_arguments(output_directory, folder_prefix_release, taxonomy):
    """
    This method convert the list of parameters into a dictionary hash that will be use by the
    corresponding pipeline.

    :param output_directory:
    :param folder_prefix_release:
    :param taxonomy:
    :return:
    """
    pipeline_arguments = {}
    if output_directory is not None:
        pipeline_arguments[EnsemblDataDownloadService._CONFIG_OUTPUT_DIRECTORY] = output_directory
    if folder_prefix_release is not None:
        pipeline_arguments[EnsemblDataDownloadService._CONFIG_KEY_FOLDER_PREFIX_RELEASE] = folder_prefix_release
    if taxonomy is not None:
        pipeline_arguments[EnsemblDataDownloadService._CONFIG_TAXONOMY] = taxonomy

    return pipeline_arguments


@click.command()
@click.option('--config_file', '-c', help='Configuration file for the ensembl data downloader pipeline',
              default='../config/config_ensembl_downloader.yaml')
@click.option('--output_directory', '-o', help='Output directory for the peptide databases', default="./database/")
@click.option('--folder_prefix_release', '-fp', help='Output folder prefix to download the data', default='release-')
@click.option('--taxonomy', '-t',
              help='Taxonomy List (comma separated) that will be use to download the data from Ensembl', default='')
@click.pass_context
def main(ctx, config_file, output_directory, folder_prefix_release, taxonomy):
    """
    The main method takes a set of parameters and create a configuration environment for the pipeline.
    :param ctx: tool context
    :param config_file: config file
    :param output_directory: output directory
    :param folder_prefix_release: folder prefix release
    :param taxonomy: taxonomies to download
    """
    if config_file is None:
        print_help(ctx, None, value=True)

    if config_file is None:
        msg = "The config file for the pipeline is missing, please provide one "
        logging.error(msg)
        raise AppConfigException(msg)

    pipeline_arguments = get_pipeline_arguments(output_directory, folder_prefix_release, taxonomy)
    ensembl_download_service = EnsemblDataDownloadService(config_file, pipeline_arguments)

    logger = ensembl_download_service.get_logger_for("Main Pipeline Ensembl Downloader")
    logger.info("Pipeline STARTING ... ")

    ensembl_download_service.download_database_by_species()

    logger.info("Pipeline Finish !!!")


if __name__ == "__main__":
    main()
