import argparse
import os
import sys

import pandas as pd

import data_processor
import modeler
import utils
import config as configuration
import json


#
# prepare --config /Users/c/Documents/umb/workspace/cs696/config.json
# extract --config /Users/c/Documents/umb/workspace/cs696/config.json --output_dir output --ids "NP_892006.3" "NP_958782.1"
#
# This is the starting point of the data preparation and processing
# for the RNA prediction application.
#
def main():
    """

    :return:
    """
    parser = argparse.ArgumentParser(prog='RNA Predictor',
                                     description='Prepare, Extract and Process RNA data',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--config',
                        required=True,
                        help='path to program config file')
    parser.add_argument('action',
                        choices=['prepare', 'extract', 'model', 'prep-single-cell',
                                 'predict-single-cell', 'calc-normalization-factor'],
                        help='action to perform\n'
                             '  prepare: process and prepare the raw proteomic and RNA data for processing\n'
                             '  extract: extract samples with non-zero rna values for a set of protein ids\n'
                             '  model: make and test prediction models based on prepared protein and rna counts\n'
                             '  calc-normalization-factor: calculate and save the single cell normalization factor\n'
                             '  prep-single-cell: prepare the single cell rna counts to for predicting protein ids\n'
                             '  predict-single-cell: predict the value of a specific protein for a specific cell\n')

    parser.add_argument('--prot_ids',
                        nargs='+',
                        help='list of protein ids to extract or predict')
    parser.add_argument('--cell_ids',
                        nargs='+',
                        help='list of cell ids to predict protein vallues')
    parser.add_argument('--output_dir',
                        nargs='?', help='directory where the results should be saved')

    args = parser.parse_args()

    # load the config json file, which contains references to the data
    if args.config is not None:
        with open(args.config) as config_file:
            config_map = json.load(config_file)
            config = configuration.Config(config_map)

    action = args.action

    #
    # calculate the normalization factor and write the value in a file to be used later
    # for prediction. This value ise used to normalize the rna values for single cells
    # to be comparable to RNA values from the proteomic data.
    #
    if action == 'calc-normalization-factor':
        with (open(config.cell_rna_total_count_file(), 'r') as cell_rna_total_count_map_file):
            cell_rna_total_count_map = json.load(cell_rna_total_count_map_file)
            normalization_factor = data_processor.find_normalization_factor(
                config.get_prepared_rna_count_file(), cell_rna_total_count_map)

            print(f'Average total normalization factor: {normalization_factor}')

            stamped_out_file = f'{config.get_normalization_fact_file()}_{utils.get_date_time_str()}'
            with open(config.get_normalization_fact_file(), 'w') as stamped_output, \
                open(stamped_out_file, 'w') as latest_output:
                stamped_output.write(str(normalization_factor))
                latest_output.write(str(normalization_factor))

    #
    # This is a preprocessing step where the relevant data are extracted from
    # the original files, aggregated into a single and patient ids are
    # reconciled from the protein and rna files.
    #
    # This step is required before any other action.
    #
    elif action == 'prepare':
        print(f'Preparing all data for processing ...')
        print(f'protemoic: {config.get_proteomic_data_dir()}')
        if not os.path.exists(config.get_proteomic_data_dir()):
            print(f'Unzipping proteomic data file {config.get_proteomic_zip_file()} ...')
            utils.unzip(config.get_proteomic_zip_file(), config.get_proteomic_data_dir())
        else:
            print(f'Found proteomic data file {config.get_proteomic_zip_file()} ...')

        prot_df, rna_df, total_cnt, na_cnt = \
            data_processor.prepare_data(raw_prot_data_dir=config.root_dir,
                                        prepared_data_dir=config.get_prepared_data_dir(),
                                        rna_counts_file=config.get_rna_counts_file(),
                                        prepared_rna_counts_file=config.get_prepared_rna_count_file(),
                                        prepared_prot_counts_file=config.get_prepared_prot_count_file())

        print(f"Across all files: total Count: {total_cnt}, NA Count: {na_cnt if na_cnt >= 0 else '<not counted>'}")

    #
    # Extract the RNA data for the specified list of proteins and save the results.
    # Only RNA samples that have non-zero values for all RNA counts are saved, that
    # is, samples that have a zero count for one or more RNA are dropped from the results.
    #
    # A list of protein ids must be provided for this action. Protein ides are stored
    # in a file called prot.csv and sample rna counts for the specified proteins are
    # stored in a file called rna.csv. both files are stored under the folder specified
    # with the --output_dir option under the program's "generated" folder. If the
    # value of --output_dir is an absolute path, the "generated" folder is ignored and
    # the results are stored under the specified folder.
    # Note that you must run the "prepare" command before running this command.
    #
    elif action == 'extract':
        # if args.prot_ids is None or args.output_dir is None:
        if args.prot_ids is None or args.output_dir is None:
            print(f'Usage: You must provide a list of protein ids and output directory with the extract command')
            parser.print_usage()
            exit(1)

        # Columns names have been normalized during the preparation step,
        # so we don't have to worry about it here

        # load the prepared data into dataframes
        df_prot_counts = pd.read_csv(config.get_prepared_prot_count_file(), sep=configuration.comma, header='infer')
        df_rna_counts = pd.read_csv(config.get_prepared_rna_count_file(), sep=configuration.comma)

        # extract the columns from the rna files that have no zeros in their values
        df_rna_counts_no_zero = df_rna_counts[df_rna_counts.iloc[:, 1] > 100]

        orig_sample_counts = df_rna_counts.shape[0]
        non_zero_sample_count = df_rna_counts_no_zero.shape[0]

        print(f'Total samples:      {orig_sample_counts}')
        print(f'Eliminated samples: {orig_sample_counts - non_zero_sample_count}')
        print(f'Extracted samples:  {non_zero_sample_count}')

        # find the common columns between the rna-non-zero-df and prot_df
        common_col = list(set(df_rna_counts_no_zero.columns).intersection(df_prot_counts.columns))
        # add back the first column name into the columns list
        filtered_prot_cols = list(common_col)
        filtered_prot_cols.insert(0, df_prot_counts.columns[0])

        filtered_rna_cols = list(common_col)
        filtered_rna_cols.insert(0, df_rna_counts.columns[0])

        filtered_prot_df = data_processor.extract_columns(df_prot_counts, filtered_prot_cols)
        filtered_rna_df = data_processor.extract_columns(df_rna_counts, filtered_rna_cols)

        # Persist the extracted data
        extract_output_dir = os.path.join(config.get_extract_data_dir(), args.output_dir)
        utils.ensure_dir_exists(extract_output_dir)
        prot_output_file = os.path.join(extract_output_dir, 'prot.csv')
        rna_output_file = os.path.join(extract_output_dir, 'rna.csv')
        print(f'writing results:\n\textracted protein counts: {prot_output_file}'
              f'\n\tExtracted RNA counts: {rna_output_file}')
        filtered_prot_df.to_csv(prot_output_file, header=True, index=False)
        filtered_rna_df.to_csv(rna_output_file, header=True, index=False)

    #
    # Create a model and test for protein prediction for a set of protein ids (specified
    # on the command line, or all protein ids in the protein count file. Results are
    # stored in a file. Models are not saved.
    # Note that running this action for all proteins will take several hours.
    #
    elif action == 'model':
        prot_ids = args.prot_ids if args.prot_ids is not None and len(args.prot_ids) > 0 else None
        results_file = modeler.test_models(config.get_prepared_prot_count_file(),
                                           config.get_prepared_rna_count_file(),
                                           config.get_model_test_results_dir(),
                                           protein_ids=prot_ids)

        print(f'Saved results in: {results_file}')

    #
    # Prepare the single cell data which can be used for predicting protein values
    # based on RNA values measured for single cells. The result of this action is
    # a json file containing a map whose keys are single cell ids and whose values
    # are maps containing RNA IDs as keys and RNA counts as values.
    # This file can then be used to create a model for a particular cell, based on
    # proteomic data, and used for predicting the count for a specific protein.
    # The model will be for a specific cell and a specific protein.
    #
    elif action == 'prep-single-cell':
        cell_rna_counts_map, non_matching_rna_names, cell_rna_total_count_map = \
            data_processor.extract_single_cell_non_zero(
                var_map_file=config.get_var_map(),
                cell_rna_id_counts_file=config.get_single_cell_rna_file(),
                normalization_method=config.get_normalization_method())
                # rna_counts_file=config.get_prepared_rna_count_file())

        # utils.log_normalize_single_cells(cell_rna_total_count_map)

        with open(config.get_single_cell_rna_counts_file(), 'w', encoding='UTF-8') as rna_counts_json, \
                open(config.get_single_cell_no_match_rna_file(), 'w', encoding='UTF-8') as non_matching_rna_json, \
                open(config.cell_rna_total_count_file(), 'w', encoding='UTF-8') as single_cell_total_count_json:
            rna_counts_json.write(json.dumps(cell_rna_counts_map, indent=2))
            non_matching_rna_json.write(json.dumps(list(non_matching_rna_names), indent=2))
            single_cell_total_count_json.write(json.dumps(cell_rna_total_count_map, indent=2))

    #
    # This command uses the RNA values for a single cell to predict the protein values for
    # the same cell. This command may be invoked for one or more proteins. For each protein
    # the program creates a model based on common RNA values with the RNA values from the
    # proteomic data set and uses that model for prediction.
    #
    elif action == 'predict-single-cell':
        if args.prot_ids is None:
            print(f'You must provide at least one protein id by using --prot_ids argument with this command')
            exit(2)

        try:
            normalization_fact_str = utils.read_line_from_file(config.get_normalization_fact_file(), 1)
        except FileNotFoundError:
            print(f'Normalization file "{config.get_normalization_fact_file()}" not found')
            sys.exit(1)

        if normalization_fact_str is None or len(normalization_fact_str) == 0:
            print(f'Normalization factor extracted from "{config.get_normalization_fact_file()}" is invalid')
            sys.exit(1)
        try:
            normalization_fact = float(normalization_fact_str)
        except ValueError:
            print(f'Normalization factor "{normalization_fact_str}" extracted from '
                  f'"{config.get_normalization_fact_file()}" is invalid')
            sys.exit(1)

        prot_ids = args.prot_ids

        cell_ids = args.cell_ids if args.cell_ids is not None and len(args.cell_ids) > 0 else None

        # this is a map of maps to hold the results of multiple protein prediction and multiple cell ids
        #   {
        #       protein_id_1 -> {
        #           cell_1: val_1
        #           ...
        #           cell_n: val_n
        #       },
        #       ...
        #       protein_id_m -> {
        #           cell_1: val_1
        #           ...
        #           cell_n: val_n
        #       }
        #
        predicted_cell_protein_count_map = {}

        # predict one protein at a time and save the results in the prediction map
        for prot_id in prot_ids:
            # the return value of this method is a map of cell id to model results,
            # which contains the predicted value for each cell-protein pair,
            single_cell_rna_count_map = data_processor.load_single_cell_counts(config.get_single_cell_rna_counts_file())

            cell_id_model_results_map = modeler.predict_prot_for_single_cell_with_df(
                prot_df=pd.read_csv(config.get_prepared_prot_count_file()),
                rna_df=pd.read_csv(config.get_prepared_rna_count_file()),
                prot_id=prot_id,
                single_cell_rna_count_map=single_cell_rna_count_map,
                normalization_fact=normalization_fact,
                cell_ids=cell_ids)

            predicted_cell_protein_count_map[prot_id] = cell_id_model_results_map

        # if an output directory is mentioned, save the results
        if args.output_dir is not None:
            output_dir = config.get_single_cell_prediction_output_dir(args.output_dir)
            utils.ensure_dir_exists(output_dir)
            output_file = os.path.join(output_dir, 'single_cell_prediction_results.json')
            with open(output_file, 'w', encoding='UTF-8') as prediction_json:
                prediction_json.write(json.dumps(predicted_cell_protein_count_map, indent=2))

            print(f'Prediction results saved in {output_file}')

        else:
            print(predicted_cell_protein_count_map)


if __name__ == "__main__":
    """
    """
    main()
