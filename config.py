import os

field_separator = '\t'
comma = ','
filtered_file_postfix = '-filtered'


#
# This class is a single place to retrieve names and locations, among other configuration
# values. There are two sets of information kept here: user defined information, such
# as the root data directory, and program defined information values such as the name
# of the generated rna counts file. The user defined information is sourced from a json
# configuration file that is specified on the command line, when the program is run.
#
# @Author: Cyrus
class Config:
    def __init__(self, config_map):
        """
        Initialize an object of the config from the specified configuration map.

        :param config_map: a dictionary object containing all user-defined values.
        """
        self.root_dir = os.path.abspath(config_map['data_root_dir'])
        self.generated_dir = os.path.abspath(config_map['generated_data_dir'])
        self.config_map = config_map

    def get_normalization_method(self):
        return self.config_map['single_cell_normalization_method']

    def get_prepared_data_dir(self):
        return os.path.join(self.generated_dir, 'prepared')

    def get_extract_data_dir(self):
        return os.path.join(self.generated_dir, 'extract')

    def get_proteomic_data_dir(self):
        return os.path.join(self.root_dir, self.config_map['proteomic_data_dir'])

    def get_proteomic_zip_file(self):
        return os.path.join(self.root_dir, self.config_map['proteomic_data_dir'] + '.zip')

    def get_rna_counts_file(self):
        return os.path.join(self.root_dir, self.config_map['tumor_rna_counts'])

    def get_var_map(self):
        return os.path.join(self.root_dir, self.config_map['var_map'])

    def get_single_cell_rna_file(self):
        return os.path.join(self.root_dir, self.config_map['single_cell_rna_file'])

    def get_prepared_prot_count_file(self):
        return os.path.join(self.get_prepared_data_dir(), 'protein_counts.csv')

    def get_prepared_rna_count_file(self):
        return os.path.join(self.get_prepared_data_dir(), 'rna_counts.csv')

    def get_filtered_dir(self):
        return os.path.join(self.generated_dir, 'filtered')

    def get_prot_file(self):
        return os.path.join(self.generated_dir, )

    def get_filtered_dir(self):
        return os.path.join(self.generated_dir, 'filtered')

    def get_single_cell_rna_counts_file(self):
        return os.path.join(self.generated_dir, 'single_cell_rna_counts.json')

    def get_single_cell_no_match_rna_file(self):
        return os.path.join(self.generated_dir, 'non_matching_rna_ids.json')

    def cell_rna_total_count_file(self):
        return os.path.join(self.generated_dir, 'single_cell_total_rna_counts.json')

    def get_model_test_results_dir(self):
        return os.path.join(self.generated_dir, 'model-test-results')

    def get_single_cell_prediction_output_dir(self, output_dir):
        return os.path.join(self.generated_dir, output_dir)

    def get_normalization_fact_file(self):
        return os.path.join(self.generated_dir, 'total_normalization_factor')
