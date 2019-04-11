"""
This module implements en Ensembl data grabber for a given Ensembl Service instance.

Some of the use cases for this module:
    1. Given a species ID, download its protein sequence data, with the option of decompressing it or not.
    2. Given a species ID, collect its GTF data, with the option of decompressing it or not.
"""

import os

# App imports
from json import loads
from requests import get
from toolbox.general import ParameterConfiguration, check_create_folders, download_file


class EnsemblDataDownloadService(ParameterConfiguration):
    """
    This Service is in charge of grabbing data (download) from Ensembl to a local repository
    """

    _CONFIG_KEY_DATA_DOWNLOADER = 'ensembl_data_downloader'
    _CONFIG_OUTPUT_DIRECTORY = 'output_directory'
    _CONFIG_KEY_ENSEMBL_FTP = 'ensembl_ftp'
    _CONFIG_ENSEMBL_API = 'ensembl_api'
    _CONFIG_ENSEMBL_API_SERVER = 'server'
    _CONFIG_ENSEMBL_API_SPECIES = 'species'
    _CONFIG_KEY_BASE_URL = 'base_url'
    _CONFIG_KEY_FOLDER_PREFIX_RELEASE = 'folder_prefix_release'
    _CONFIG_KEY_FOLDER_NAME_FASTA = 'folder_name_fasta'
    _CONFIG_KEY_FOLDER_NAME_PROTEIN_SEQUENCES = 'folder_name_protein_sequences'
    _CONFIG_KEY_FOLDER_NAME_GTF = 'folder_name_gtf'
    _CONFIG_KEY_REWRITE_LOCAL_PATH_ENSEMBL_REPO = 'rewrite_local_path_ensembl_repo'
    _CONFIG_KEY_ENSEMBL_FILE_NAMES = 'ensembl_file_names'
    _CONFIG_KEY_PROTEIN_SEQUENCE_FILE = 'protein_sequence_file'
    _CONFIG_KEY_FILE_TYPE = 'file_type'
    _CONFIG_KEY_FILE_SUFFIXES = 'file_suffixes'
    _CONFIG_KEY_FILE_EXTENSION = 'file_extension'
    _CONFIG_KEY_GTF_FILE = 'gtf_file'
    _CONFIG_REST_API_TAXON_ID = 'taxon_id'
    _CONFIG_TAXONOMY = 'taxonomy'

    def __init__(self, config_file, pipeline_arguments):
        """
        Init the class with the specific parameters.
        :param params:
        """
        super(EnsemblDataDownloadService, self).__init__(self._CONFIG_KEY_DATA_DOWNLOADER, config_file,
                                                         pipeline_arguments)

        self._ensembl_species = []
        if self._CONFIG_OUTPUT_DIRECTORY in self.get_pipeline_parameters():
            self._local_path_ensembl = self.get_pipeline_parameters()[self._CONFIG_OUTPUT_DIRECTORY]
        else:
            self._local_path_ensembl = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                self._CONFIG_OUTPUT_DIRECTORY]

        self.prepare_local_ensembl_repository()

    def get_local_path_root_ensembl_repo(self):
        return self._local_path_ensembl

    def prepare_local_ensembl_repository(self):
        self.get_logger().debug(
            "Preparing local Ensembl repository, root folder - '{}'".format(self.get_local_path_root_ensembl_repo()))
        check_create_folders([self.get_local_path_root_ensembl_repo()])
        self.get_logger().debug(
            "Local path for Ensembl Release - '{}'".format(self.get_local_path_root_ensembl_repo()))

    def get_species_from_rest(self):
        """
        Get the list of species from ENSEMBL rest API.
        :return:
        """
        server = self.get_default_parameters()[self._CONFIG_KEY_DATA_DOWNLOADER][self._CONFIG_ENSEMBL_API][
            self._CONFIG_ENSEMBL_API_SERVER]
        endpoint = self.get_default_parameters()[self._CONFIG_KEY_DATA_DOWNLOADER][self._CONFIG_ENSEMBL_API][
            self._CONFIG_ENSEMBL_API_SPECIES]
        species_info = loads(get(server + endpoint, headers={"Content-Type": 'application/json'}).text)
        self._ensembl_species = species_info['species']
        return self._ensembl_species

    def download_database_by_species(self):
        """
        This method takes a list of Taxonomies from the commandline parameters and download the Protein fasta files
        and the gtf files.
        :return:
        """
        self.get_species_from_rest()
        species_parameters = self.get_pipeline_parameters()[self._CONFIG_TAXONOMY]
        species_list = species_parameters.split(",")
        total_files = []
        if species_list is None or len(species_list) == 0 or len(species_parameters) == 0:
            for species in self._ensembl_species:
                self.get_logger().debug(
                    "Downloading the data for the specie -- " + species[self._CONFIG_REST_API_TAXON_ID])
                files = self.get_pep_files(species)
                gtf_files = self.gt_gtf_files(species)
                files.extend(gtf_files)
                total_files.extend(files)
                self.get_logger().debug("Files downloaded -- " + ",".join(files))
                total_files.extend(files)
        else:
            for species_id in species_list:
                for species in self._ensembl_species:
                    if species_id == species[self._CONFIG_REST_API_TAXON_ID]:
                        self.get_logger().debug(
                            "Downloading the data for the specie -- " + species[self._CONFIG_REST_API_TAXON_ID])
                        files = self.get_pep_files(species)
                        gtf_files = self.gt_gtf_files(species)
                        files.extend(gtf_files)
                        total_files.extend(files)
                        self.get_logger().debug("Files downloaded -- " + ",".join(files))
                        total_files.extend(files)

        return total_files

    def get_pep_files(self, species: dict) -> list:
        """
        Get the peptide files for an specific species object.
        :return: List of files names.
        """
        files = []
        try:
            file_name = '{}.{}.pep.all.fa.gz'.format(species['name'][0].upper() + species['name'][1:],
                                                     species['assembly'])
            file_url = '{}/release-{}/fasta/{}/pep/{}'.format(
                self.get_default_parameters()[self._CONFIG_KEY_DATA_DOWNLOADER][self._CONFIG_KEY_ENSEMBL_FTP][
                    self._CONFIG_KEY_BASE_URL],
                species['release'], species['name'], file_name)
            files.append(download_file(file_url, self.get_local_path_root_ensembl_repo() + '/' + file_name))
        except KeyError:
            print("No valid info is available species: ", species)

        return files

    def gt_gtf_files(self, species: dict) -> list:
        """
        This method retrieve the gtf files for an specific specie object
        :param species:
        :return:
        """
        """
          Generate GTF file name from the species info and download the GTF file
          """
        files = []
        try:
            file_name = '{}.{}.{}.gtf.gz'.format(species['name'][0].upper() + species['name'][1:], species['assembly'],
                                                 species['release'], )
            file_url = '{}/release-{}/gtf/{}/{}'.format(
                self.get_default_parameters()[self._CONFIG_KEY_DATA_DOWNLOADER][self._CONFIG_KEY_ENSEMBL_FTP][
                    self._CONFIG_KEY_BASE_URL], species['release'], species['name'], file_name)
            files.append(download_file(file_url, self.get_local_path_root_ensembl_repo() + '/' + file_name))
        except KeyError:
            print("No valid info is available species: ", species)

        return files


if __name__ == '__main__':
    print("ERROR: This script is part of a pipeline collection and it is not meant to be run in stand alone mode")
