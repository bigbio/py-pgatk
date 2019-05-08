from Bio import SeqIO

from toolbox.general import ParameterConfiguration


class EnsemblDataService(ParameterConfiguration):

    CONFIG_KEY_DATA = "enmsembl_translation"
    CONFIG_TRANSLATION_TABLE = "translation_table"

    def __init__(self, config_file, pipeline_arguments):
        """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
        super(EnsemblDataService, self).__init__(self.CONFIG_KEY_DATA, config_file,
                                                 pipeline_arguments)

        if self.CONFIG_TRANSLATION_TABLE in self.get_pipeline_parameters():
            self._translation_table = self.get_pipeline_parameters()[self.CONFIG_TRANSLATION_TABLE]
        else:
            self._translation_table = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_TRANSLATION_TABLE]

    def three_frame_translation(self, input_file, output_file):

        input_handle = open(input_file, 'r')
        output_handle = open(output_file, 'rw')

        for record in SeqIO.parse(input_handle, 'fasta'):
            seq = record.seq
            RF1 = seq.translate(table=self._translation_table)
            RF2 = seq[1::].translate(table=self._translation_table)
            RF3 = seq[2::].translate(table=self._translation_table)

            if record.id == "":
                print("skip entries without id", record.description)
                continue
            output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF1', RF1))
            output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF2', RF2))
            output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF3', RF3))


if __name__ == '__main__':
    print("ERROR: This script is part of a pipeline collection and it is not meant to be run in stand alone mode")
