import logging

import click

from pypgatk.proteomics.db.decoy_pyrat import ProteinDBService
from pypgatk.toolbox.exceptions import AppConfigException


@click.command('proteindb-decoy', short_help='Create decoy protein sequences. Each protein '
                                             'is reversed and the cleavage sites switched with preceding amino acid. Peptides are checked for existence in target sequences if found'
                                             'the tool will attempt to shuffle them. James.Wright@sanger.ac.uk 2015')
@click.option('--config_file', '-c', help='Configuration file for the protein database decoy generation',
              default='config/protein_decoy.yaml')
@click.option('--output', '-o', help='Output file for decoy database',
              default="protein-decoy.fa")
@click.option('--input','-i',  metavar='*.fasta|*.fa', help='FASTA file of target proteins sequences for which to'
                                                            ' create decoys')
@click.option('--cleavage_sites', '-c', default='KR',
                    help='A list of amino acids at which to cleave during digestion. Default = KR')
@click.option('--anti_cleavage_sites', '-a', dest='noc', default='',
                    help='A list of amino acids at which not to cleave if following cleavage site ie. Proline. Default = none')
@click.option('--cleavage_position', '-p', dest='cpos', default='c', choices=['c', 'n'],
                    help='Set cleavage to be c or n terminal of specified cleavage sites. Default = c')
@click.option('--min_peptide_length', '-l', dest='minlen', default=5, type=int,
                    help='Set minimum length of peptides to compare between target and decoy. Default = 5')
@click.option('--max_iterations', '-n', dest='maxit', default=100, type=int,
                    help='Set maximum number of times to shuffle a peptide to make it non-target before failing. Default=100')
@click.option('--do_not_shuffle', '-x', dest='noshuf', default=False, action='store_true',
                    help='Turn OFF shuffling of decoy peptides that are in the target database. Default=false')
@click.option('--do_not_switch', '-s', dest='noswitch', default=False, action='store_true',
                    help='Turn OFF switching of cleavage site with preceding amino acid. Default=false')
@click.option('--decoy_prefix', '-d', dest='dprefix', default='XXX',
                    help='Set accession prefix for decoy proteins in output. Default=XXX')
@click.option('--temp_file', '-t', dest='tout', default='tmp.fa',
                    help='Set temporary file to write decoys prior to shuffling. Default=tmp.fa')
@click.option('--no_isobaric', '-b', dest='iso', default=False, action='store_true',
                    help='Do not make decoy peptides isobaric. Default=false')
@click.option('--memory_save', '-m', dest='mem', default=False, action='store_true',
                    help='Slower but uses less memory (does not store decoy peptide list). Default=false')
@click.pass_context
def decoy_database(ctx, config_file, output, input, cleavage_sites, anti_cleavages_sites, cleavage_position, min_peptide_length,
                   max_interactions, do_not_shuffle, do_not_switch, decoy_prefix, temp_file, no_isobaric, memory_save):

    if config_file is None:
        msg = "The config file for the pipeline is missing, please provide one "
        logging.error(msg)
        raise AppConfigException(msg)

    pipeline_arguments = {}
    if output is not None:
        pipeline_arguments[ProteinDBService.CONFIG_OUTPUT_DIRECTORY] = output
    if input is not None:
        pipeline_arguments[ProteinDBService.INPUT_FILE] = input

    cbioportal_downloader_service = CbioPortalDownloadService(config_file, pipeline_arguments)

    if list_studies:
        list_studies = cbioportal_downloader_service.get_cancer_studies()
        print(list_studies)

    if download_study is not None:
        cbioportal_downloader_service.download_study(download_study)