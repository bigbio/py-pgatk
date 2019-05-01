import base64
import json

import requests
from toolbox.general import ParameterConfiguration, check_create_folders


class CosmicDownloadService(ParameterConfiguration):
    _CONFIG_KEY_DATA_DOWNLOADER = 'cosmic_data_downloader'
    _CONFIG_OUTPUT_DIRECTORY = 'output_directory'
    _CONFIG_COSMIC_SERVER = 'cosmic_server'
    _CONFIG_COSMIC_FTP_URL = 'cosmic_ftp'
    _CONFIG_COSMIC_FTP_USER = "cosmic_user"
    _CONFIG_COSMIC_FTP_PASSWORD = "cosmic_password"
    _CONFIG_COSMIC_MUTATIONS_URL = "mutations_url"
    _CONFIG_COSMIC_MUTATIONS_FILE = "mutations_file"

    def __init__(self, config_file, pipeline_arguments):
        """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
        super(CosmicDownloadService, self).__init__(self._CONFIG_KEY_DATA_DOWNLOADER, config_file, pipeline_arguments)

        if self._CONFIG_OUTPUT_DIRECTORY in self.get_pipeline_parameters():
            self._local_path_cosmic = self.get_pipeline_parameters()[self._CONFIG_OUTPUT_DIRECTORY]
        else:
            self._local_path_cosmic = self.get_default_parameters()[self._CONFIG_KEY_DATA_DOWNLOADER][
                self._CONFIG_OUTPUT_DIRECTORY]

        self._cosmic_token = base64.b64encode(
            f'{self.get_pipeline_parameters()[self._CONFIG_COSMIC_FTP_USER]}:{self.get_pipeline_parameters()[self._CONFIG_COSMIC_FTP_PASSWORD]}'.encode()).decode(
            'utf-8')

        self.prepare_local_cosmic_repository()

    def prepare_local_cosmic_repository(self):
        self.get_logger().debug("Preparing local cbioportal repository, root folder - '{}'".format(
            self.get_local_path_root_cosmic_repo()))
        check_create_folders([self.get_local_path_root_cosmic_repo()])
        self.get_logger().debug(
            "Local path for cbioportal Release - '{}'".format(self.get_local_path_root_cosmic_repo()))

    def get_local_path_root_cosmic_repo(self):
        return self._local_path_cosmic

    def download_mutation_file(self):
        """
        This function will download the mutations file from Cosmic Database.
        :return: None
        """

        output_file = f'{self.get_local_path_root_cosmic_repo()}/{self.get_default_parameters()[self._CONFIG_KEY_DATA_DOWNLOADER][self._CONFIG_COSMIC_SERVER][self._CONFIG_COSMIC_MUTATIONS_FILE]}'

        server = self.get_default_parameters()[self._CONFIG_KEY_DATA_DOWNLOADER][self._CONFIG_COSMIC_SERVER][
            self._CONFIG_COSMIC_FTP_URL]
        cosmic_version = self.get_default_parameters()[self._CONFIG_KEY_DATA_DOWNLOADER][self._CONFIG_COSMIC_SERVER][
            self._CONFIG_COSMIC_MUTATIONS_URL]
        mutation_file = self.get_default_parameters()[self._CONFIG_KEY_DATA_DOWNLOADER][self._CONFIG_COSMIC_SERVER][
            self._CONFIG_COSMIC_MUTATIONS_FILE]
        url = f'{server}/{cosmic_version}/{mutation_file}'

        token = f'Basic {self._cosmic_token}'
        response = requests.get(url, stream=True, headers={'Authorization': token})

        if response.status_code == 200:
            url = json.loads(response.text)['url']
            msg = "Downloading file from url '{}'".format(url)
            self.get_logger().debug(msg)

            response = requests.get(url, stream=True)

            if response.status_code == 200:
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                    msg = "Download Finish for file '{}'".format(output_file)
                    self.get_logger().debug(msg)
