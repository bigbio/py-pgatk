import csv

from requests import get

from pypgatk.toolbox.exceptions import AppException
from pypgatk.toolbox.general import ParameterConfiguration, check_create_folders, download_file


class CbioPortalDownloadService(ParameterConfiguration):
    CONFIG_KEY_DATA_DOWNLOADER = 'cbioportal_data_downloader'
    CONFIG_KEY_CBIOPORTAL_DOWNLOAD_URL = 'cbioportal_download_url'
    CONFIG_OUTPUT_DIRECTORY = 'output_directory'
    CONFIG_CBIOPORTAL_API = 'cbioportal_api'
    CONFIG_CBIOPORTAL_API_SERVER = 'base_url'
    CONFIG_CBIOPORTAL_API_CANCER_STUDIES = "cancer_studies"
    CONFIG_LIST_STUDIES = "list_studies"

    def __init__(self, config_file, pipeline_arguments):
        """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
        super(CbioPortalDownloadService, self).__init__(self.CONFIG_KEY_DATA_DOWNLOADER, config_file,
                                                        pipeline_arguments)

        self.cbioportal_studies = []
        if self.CONFIG_OUTPUT_DIRECTORY in self.get_pipeline_parameters():
            self._local_path_cbioportal = self.get_pipeline_parameters()[self.CONFIG_OUTPUT_DIRECTORY]
        else:
            self._local_path_cbioportal = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                self.CONFIG_OUTPUT_DIRECTORY]

        self.prepare_local_cbioportal_repository()

    def prepare_local_cbioportal_repository(self):
        self.get_logger().debug("Preparing local cbioportal repository, root folder - '{}'".format(
            self.get_local_path_root_cbioportal_repo()))
        check_create_folders([self.get_local_path_root_cbioportal_repo()])
        self.get_logger().debug(
            "Local path for cbioportal Release - '{}'".format(self.get_local_path_root_cbioportal_repo()))

    def get_local_path_root_cbioportal_repo(self):
        return self._local_path_cbioportal

    def get_cancer_studies(self):
        """
        This method will print the list of all cancer studies for the user.
        :return:
        """
        server = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_CBIOPORTAL_API][
            self.CONFIG_CBIOPORTAL_API_SERVER]
        endpoint = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_CBIOPORTAL_API][
            self.CONFIG_CBIOPORTAL_API_CANCER_STUDIES]
        self._cbioportal_studies = get(server + "?" + endpoint).text
        return self._cbioportal_studies

    def download_study(self, download_study):
        """
        This function will download an study from cBioPortal using the study ID
        :param download_study: Study to be download, if the study is empty or None, all the studies will be
        download.
        :return: None
        """

        if self._cbioportal_studies is None or len(self._cbioportal_studies):
            self.get_cancer_studies()

        if 'all' not in download_study:
            if not self.check_study_identifier(download_study):
                msg = "The following study accession '{}' is not present in cBioPortal Studies".format(download_study)
                self.get_logger().debug(msg)
                raise AppException(msg)
            else:
                self.download_one_study(download_study)
        else:
            csv_reader = csv.reader(self._cbioportal_studies.splitlines(), delimiter="\t")
            line_count = 0
            for row in csv_reader:
                if line_count != 0:
                    self.download_one_study(row[0])
                line_count = line_count + 1

    def download_one_study(self, download_study):
        file_name = '{}.tar.gz'.format(download_study)
        file_url = '{}/{}'.format(
            self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_KEY_CBIOPORTAL_DOWNLOAD_URL],
            file_name)
        file_name = download_file(file_url, self.get_local_path_root_cbioportal_repo() + '/' + file_name)
        msg = "The following study '{}' has been downloaded. ".format(file_name)
        self.get_logger().debug(msg)
        return file_name

    def check_study_identifier(self, download_study):
        return download_study in self._cbioportal_studies
