from toolbox.general import ParameterConfiguration, check_create_folders


class CbioPortalDownloadService(ParameterConfiguration):


    _CONFIG_KEY_DATA_DOWNLOADER = 'cbioportal_data_downloader'
    _CONFIG_OUTPUT_DIRECTORY = 'output_directory'
    _CONFIG_CBIOPORTAL_API = 'cbioportal_api'
    _CONFIG_CBIOPORTAL_API_SERVER = 'base_url'

    def __init__(self, config_file, pipeline_arguments):
        """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
        super(CbioPortalDownloadService, self).__init__(self._CONFIG_KEY_DATA_DOWNLOADER, config_file,
                                                        pipeline_arguments)

        self.cbioportal_studies = []
        if self._CONFIG_OUTPUT_DIRECTORY in self.get_pipeline_parameters():
            self._local_path_cbioportal = self.get_pipeline_parameters()[self._CONFIG_OUTPUT_DIRECTORY]
        else:
            self._local_path_cbioportal = self.get_default_parameters()[self._CONFIG_KEY_DATA_DOWNLOADER][
                self._CONFIG_OUTPUT_DIRECTORY]

        self.prepare_local_cbioportal_repository()

    def prepare_local_cbioportal_repository(self):
        self.get_logger().debug("Preparing local cbioportal repository, root folder - '{}'".format(
            self.get_local_path_root_cbioportal_repo()))
        check_create_folders([self.get_local_path_root_cbioportal_repo()])
        self.get_logger().debug("Local path for cbioportal Release - '{}'".format(self.get_local_path_root_cbioportal_repo()))

    def get_local_path_root_cbioportal_repo(self):
        return self._local_path_cbioportal
