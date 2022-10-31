import logging

import click

from pypgatk.toolbox.general import read_yaml_from_file
from pypgatk.commands.utils import print_help
from pypgatk.proteogenomics.blast_get_position import BlastGetPositionService

log = logging.getLogger(__name__)

@click.command('blast_get_position', short_help='Blast peptide and refence protein database to find variation sites.')
@click.option('-c', '--config_file', help='Configuration file for the fdr peptides pipeline.')
@click.option('-i', '--input_psm_to_blast', help='The file name of the input PSM table to blast.')
@click.option('-o', '--output_psm', help='The file name of the output PSM table.')
@click.option('-r', '--input_reference_database', help='The file name of the refence protein database to blast. The reference database includes Uniprot Proteomes with isoforms, ENSEMBL, RefSeq, etc.')
@click.option('-p', '--canonical_peptide_prefix', help='The prefix of reference database. Uniprot is prefixed with "sp",ENSEMBL is prefixed with "ENSP", RefSeq is prefixed with "NP". Default is "sp,NP,ENSP".')

@click.pass_context
def blast_get_position(ctx, config_file, input_psm_to_blast, output_psm, input_reference_database, canonical_peptide_prefix):
    
    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    if input_psm_to_blast is None or input_reference_database is None or output_psm is None:
        print_help()
    
    pipeline_arguments = {}

    if input_reference_database is not None:
        pipeline_arguments[BlastGetPositionService.CONFIG_INPUT_REFERENCE_DATABASE] = input_reference_database
    if canonical_peptide_prefix is not None:
        pipeline_arguments[BlastGetPositionService.CONFIG_CANONICAL_PEPTIDE_PREFIX] = canonical_peptide_prefix

    blast_get_position_service = BlastGetPositionService(config_data, pipeline_arguments)
    blast_get_position_service.blast(input_psm_to_blast, output_psm)