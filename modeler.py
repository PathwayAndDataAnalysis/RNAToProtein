import datetime
import os
import time
import pandas as pd
import pickle

import data_processor
import numpy as np

from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler

import df_utils
import utils

#
# This file contains functions and classes that are used for building regression models
# and predicting protein values for single cells based on their RNA count measurement.
#
# @author: Cyrus
#


#
# This is a container class to hold the model, its properties and the prediction results.
#
class ModelResults:
    def __init__(self, prot_id):
        self.lr_model = None
        self.rna_list = None
        self.cell_id = None
        self.prot_id = prot_id
        self.selected_sample_cnt = 0
        self.train_cnt = 0
        self.test_cnt = 0
        self.test_mse = 0
        self.test_r2_score = 0
        self.tr_mse = 0
        self.tr_r2_score = 0
        self.total_count = 0
        self.predicted_value = 0

    def to_csv(self):
        return f'{self.prot_id},{self.selected_sample_cnt},{self.train_cnt},{self.test_cnt},' \
               f'{self.test_mse},{self.test_r2_score},{self.tr_mse},{self.tr_r2_score}'

    def __str__(self):
        mod_mse = '{:.2f}'.format(self.test_mse)
        mod_r2_score = '{:.2f}'.format(self.test_r2_score)
        mod_mse_tr = '{:.2f}'.format(self.tr_mse)
        mod_r2_score_tr = '{:.2f}'.format(self.tr_r2_score)
        return f'Protein id: {self.prot_id}\n\ttrain count: {self.train_cnt}\n\ttest count: {self.test_cnt}\n\t' \
               f'mse: {mod_mse}\n\tr2_score: {mod_r2_score}\n\t' \
               f'mse_tr: {mod_mse_tr}\n\tr2_score_tr: {mod_r2_score_tr}'

    #
    # Save the model on disk, in the specified file, to be loaded and used later
    #
    def save_model(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self.lr_model, f)

    #
    # Load a model save in the specified path and set it on this object
    #
    def load_model(self, path):
        with open(path, 'rb') as f:
            self.lr_model = pickle.load(f)

    def to_dict(self):
        return self.__dict__


#
# This is to hold the results of a single prediction, which is for a (protein id, single cell id) pair.
# It holds the predicted value as well as some information about the model and its performance.
#
class CellPrediction:
    def __init__(self, model_results):
        self.cell_id = model_results.cell_id
        self.train_cnt = model_results.train_cnt
        self.test_cnt = model_results.test_cnt
        self.test_mse = model_results.test_mse
        self.test_r2_score = model_results.test_r2_score
        self.tr_mse = model_results.tr_mse
        self.tr_r2_score = model_results.tr_r2_score
        self.predicted_value = model_results.predicted_value


def csv_result_header():
    """
    Prepare a string containing the header for the model prediction result
    :return: a csv string for the model prediction results
    """
    return 'prot_id,selected_sample_cnt,train_cnt,self.test_cnt,test_mse,test_r2_score,train_mse,train_r2_score'


def prepare_data(rna_df, prot_df, protein_id, rna_ids, log_normalize, test_portion=0.3):
    """
    Prepare training and test data sets for protein prediction model based on RNA counts. The features
    are RNA counts and dependent variables (to be predicted)
    :param rna_df: dataframe for the RNA counts, columns are RNA values for each RNA
    :param prot_df: dataframe for the protein counts, columns are protein values for each RNA
    :param protein_id: the id of the protein for which a model is to be created
    :param rna_ids: list of RNA ids to include in the sample, None indicates that all RNAs should be used
    :param log_normalize: indicates whether the RNA values should be log normalized
    :param test_portion: the portion of the data that should be used for testing. 1 minus test_portion is used
    for training and must be between 0.1 and 0.6 inclusively
    :return: a seven-tuple containing feature values, dependent values, feature values for training, dependent values
    for training, feature values for testing, dependent values for testing, an ordered list of rna ids used in the daya
    """
    selected_rna_df = data_processor.extract_rna_samples(rna_df, rna_ids)
    selected_prot_df = data_processor.extract_protein_samples(prot_df, [protein_id])
    # drop all columns with NA value
    selected_prot_df = selected_prot_df.dropna(axis=1)
    selected_prot_df_c, selected_rna_df_c = \
        data_processor.unify_prot_and_rna_data(selected_prot_df, selected_rna_df, normalize_col_names=False)

    if log_normalize:
        selected_rna_df_c = df_utils.log_normalize_df(selected_rna_df_c, skip_cols_count=1, replace_inf_with=0)

    # drop the protein_id and RNA name from the dataframe
    y = selected_prot_df_c.drop(['protein_id'], axis=1).to_numpy().transpose()
    x = selected_rna_df_c.drop(['Name'], axis=1).to_numpy().transpose()
    ordered_rna_ids = selected_rna_df_c['Name']

    std_scaler = StandardScaler()
    x = std_scaler.fit_transform(x, y)

    samples_count = y.shape[0]

    # put some guard rails for the test portion
    if type(test_portion) is not float or test_portion < 0.1 or test_portion > 0.6:
        raise ValueError('The value of test portion argument must be a floating point between 0.1 and 0.6')

    test_count = int(samples_count * test_portion)
    # train_count = samples_count - test_count

    x_train = x[0: -test_count]
    x_test = x[-test_count:]
    y_train = y[0: -test_count]
    y_test = y[-test_count:]

    return x, y, x_train, y_train, x_test, y_test, ordered_rna_ids


def make_linear_regression(rna_df, prot_df, prot_id, rna_ids=None):
    """
    Make a linear regression model for a particular protein id using a subset of or all RNA ids.

    :param rna_df: dataframe containing the rna counts
    :param prot_df: dataframe containing the protein counts
    :param prot_id: the target protein id
    :param rna_ids: a list of RNA ids to use for the model. In None, all RNAs are used
    :return: A ModelResults object which contains the regression model and other related information
    about the model.
    """

    X, y, x_train, y_train, x_test, y_test, selected_rna_ids = (
        prepare_data(rna_df, prot_df, prot_id, rna_ids=rna_ids, log_normalize=True))

    results = ModelResults(prot_id)
    # if we're making models with all the RNAs for multiple cells, saving the list of
    # RNA ids takes up too much memory and may not be necessary
    results.rna_list = selected_rna_ids
    results.selected_sample_cnt = y.shape[0]
    results.train_cnt = y_train.shape[0]
    results.test_cnt = y_test.shape[0]

    # Create linear regression object
    regr_model = linear_model.LinearRegression()

    # Train the model using the training sets
    regr_model.fit(x_train, y_train)
    results.lr_model = regr_model

    # Make predictions using the testing set
    y_pred = regr_model.predict(x_test)
    y_pred_tr = regr_model.predict(x_train)

    # The coefficients
    # print("CoeffiCoefficients: \n", regr.coef_)

    results.test_mse = mean_squared_error(y_test, y_pred)
    results.test_r2_score = r2_score(y_test, y_pred)

    results.tr_mse = mean_squared_error(y_train, y_pred_tr)
    results.tr_r2_score = r2_score(y_train, y_pred_tr)
    return results


def extract_rna_counts(cell_rna_count_map, rna_list):
    """
    Extract from the map value for the RNAs that are contained in the specified list.
    :param cell_rna_count_map: map of rna ids and count
    :param rna_list: list of rna ids whose values to be extracted from the map
    :return: list of values for the rna ids in the specified list in the same order
    """
    counts = []
    for rna in rna_list:
        counts.append(cell_rna_count_map.get(rna))
    return counts


def adjust_single_cell_rna_count(cell_rna_count_map, normalization_fact):
    """
    multiple every value of for every key in the cell_rna_count_map by the specified factor.

    :param cell_rna_count_map: single cell rna count map
    :param normalization_fact: factor by which values are to be multiplied
    :return: void
    """
    for key in cell_rna_count_map.keys():
        cell_rna_count_map[key] = cell_rna_count_map[key] * normalization_fact


def predict_prot_for_single_cell(prot_count_file, rna_count_file, prot_id, single_cell_count_file,
                                 normalization_fact, cell_ids=None):
    """
`   Predict the value of the protein with the specified id for each one of the cells specified
    in call_ids. If cell_ids is None, the protein is predicted for all cells that are found
    in the single_cell_count_file. A linear regression model is build and trained for each
    combination of (cell_id, protein_id) and used to make the prediction.

    :param prot_count_file: path to the protein counts file
    :param rna_count_file: path to RNA counts file
    :param prot_id: protein id whose value to be predicted
    :param single_cell_count_file: path to single cell RNA count file
    :param normalization_fact: normalization factor to adjust (multiply) the single cell values by
    :param cell_ids: list of cell ids for which the value of the specified protein it to be predicted
    :return: a map of cell id to model results, which contains the predicted value for each cell-protein pair,
    as well as additional information about the model and its performance
    """
    prot_df = pd.read_csv(prot_count_file)
    rna_df = pd.read_csv(rna_count_file)
    single_cell_rna_count_map = data_processor.load_single_cell_counts(single_cell_count_file)

    return predict_prot_for_single_cell_with_df(prot_df,
                                                rna_df,
                                                prot_id,
                                                single_cell_rna_count_map,
                                                normalization_fact,
                                                cell_ids)


def predict_prot_for_single_cell_with_df(prot_df, rna_df, prot_id, single_cell_rna_count_map,
                                         normalization_fact, cell_ids=None):
    """
`   Predict the value of the protein with the specified id for each one of the cells specified
    in call_ids. If cell_ids is None, the protein is predicted for all cells that are found
    in the single_cell_count_file. A linear regression model is build and trained for each
    combination of (cell_id, protein_id) and used to make the prediction.

    :param prot_df: dataframe containing protein counts
    :param rna_df: dataframe containing rna counts
    :param prot_id: protein id whose value to be predicted
    :param single_cell_rna_count_map: a map containing counts for each rna for each cell id
    :param normalization_fact: normalization factor to adjust (multiply) the single cell values by
    :param cell_ids: list of cell ids for which the value of the specified protein it to be predicted
    :return: a map of cell id to model results, which contains the predicted value for each cell-protein pair,
    as well as additional information about the model and its performance
    """

    predictions = []
    # extract the RNA ids for the specified cell_id
    # single_cell_rna_count_map = data_processor.load_single_cell_counts(single_cell_count_file)

    cell_ids = cell_ids if cell_ids is not None and len(cell_ids) > 0 else single_cell_rna_count_map.keys()

    # for each cell_id in the list of cells_ids, make a model and run the prediction.
    # a model is built and trained for each pair of (cell_id, protein_id) because different
    # cells have different missing values for rna counts.
    for cell_id in cell_ids:
        if cell_id not in single_cell_rna_count_map:
            print(f'Cell {cell_id} not found in the counts map')
            continue
        # extract the rna counts for the cell id
        single_cell_rna_counts = single_cell_rna_count_map[cell_id]
        rna_ids = single_cell_rna_counts.keys()
        model_results = make_linear_regression(rna_df=rna_df, prot_df=prot_df, prot_id=prot_id, rna_ids=rna_ids)
        cell_rna_counts = extract_rna_counts(single_cell_rna_counts, model_results.rna_list)

        # multiply the rna counts for the single cell by the normalization factor
        cell_rna_counts = [count * normalization_fact for count in cell_rna_counts]
        values = np.array(cell_rna_counts).reshape(1, -1)

        predicted_value = model_results.lr_model.predict(values)
        model_results.cell_id = cell_id
        model_results.predicted_value = predicted_value[0][0]
        # add the model results, which contains additional information about the model
        # and its performance in the prediction list as a map
        cp = CellPrediction(model_results)
        predictions.append(cp.__dict__)
    return predictions


def test_models(get_prepared_prot_count_file, get_prepared_rna_count_file, results_dir, protein_ids=None):
    """
    This method takes the path to the prepared protein and rna count files, a path output directory
    and an optional list of protein ids and does the following:
        - prepares a model for each protein in the protein_ids list if it's not None. If it's None
        - it is initialized to all protein ids from the prepared protein count file.
        - executes the model and measures mean squared error and r2 score for each model.
        - save the results in a file in the specified results output directory for each protein
          model. It does not save the model itself.
        - keep track of the model with the minimum mse and maximum r2 scroe.
    :param get_prepared_prot_count_file: path to the protein count file
    :param get_prepared_rna_count_file: path to the ran count file
    :param results_dir: path to the directory where results should be stored
    :param protein_ids: a list of ids to model, if None, all proteins are modeled and tested
    :return: path to the results file
    """
    now = datetime.datetime.now()
    date_str = f'{now.year}_{now.month}_{now.day}_{now.hour}_{now.minute}_{now.second}'
    utils.ensure_dir_exists(results_dir)
    results_file = os.path.join(results_dir, f'lr_results_partial_{date_str}')

    with open (results_file, 'w') as r_file:
        r_file.write(f'{csv_result_header()}\n')

        prot_df = df_utils.read_df(get_prepared_prot_count_file)
        rna_df = df_utils.read_df(get_prepared_rna_count_file)

        # if no protein ids are specified, run the program for all proteins from
        # the protein count file
        if protein_ids is None:
            protein_ids = prot_df['protein_id'].tolist()

        # keep track of some stats about the results
        min_mse = 1000000
        max_r2_score = -1000000
        min_mse_prot_id = None;
        max_r2_score_prot_id = None

        protein_count = len(protein_ids)
        print(f'Total proteins: {protein_count}')
        init_start_time = time.time()
        i = 0

        for prot_id in protein_ids:
            i += 1
            start_time = time.time()
            result = make_linear_regression(rna_df, prot_df, prot_id=prot_id, rna_ids=None)
            r_file.write(f'{result.to_csv()}\n')

            if result.test_r2_score > max_r2_score:
                max_r2_score = result.test_r2_score
                max_r2_score_prot_id = result.prot_id

            if result.test_mse < min_mse:
                min_mse = result.test_mse
                min_mse_prot_id = result.prot_id

            exec_time = '{:.2f}'.format(time.time() - start_time)
            print(f'iteration {i} of {protein_count}: {exec_time} '
                  f'secs, total run time: {int(time.time() - init_start_time)}')
            print(result)
            print(f'Min mse observed so far      = {min_mse}  for protein {min_mse_prot_id}')
            print(f'Max r2_score observed so far = {max_r2_score} for protein {max_r2_score_prot_id}')
            print('___________________________________________________________________________________________________')

            estimated_remaining_secs = int(((time.time() - init_start_time) / i) * protein_count)

            print(f'Estimated Remaining Time: {utils.convert_time_in_hours_mins(estimated_remaining_secs)}')

        r_file.close()
        print(f'Ran {len(protein_ids)} models in {int(time.time() - init_start_time)} secs')

    return results_file
