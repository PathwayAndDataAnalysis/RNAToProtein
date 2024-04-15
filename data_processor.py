import os
import json
import pandas as pd

import df_utils
import utils
import statistics as stats

from df_utils import COMMA
from df_utils import FIELD_SEPARATOR
FILTERED_FILE_POSTFIX = '-filtered'


#
# This file contains functions that are used for reading and processing data, related
# to protein and rna counts.
#
# @author: Cyrus
#

def extract_rna_samples(rna_df, rna_names=None):
    """
    Extract from the specified dataframe RNA records whose name match any of the
    specified names. If the RNA name list is None, the original dataframe is returned.

    :param rna_df: dataframe containing the records
    :param rna_names: list of names of the RNAs to match
    :return: a dataframe containing records whose name is in the specified list
    """
    return df_utils.filter_dataframe(rna_df, col_name='Name', values=rna_names)


def extract_protein_samples(prot_df, prot_ids=None):
    """
    Extract from the specified dataframe protein records whose name match any of the
    specified names. If the protein ids list is None, the original dataframe is returned.

    :param prot_df: dataframe containing the records
    :param prot_ids: list of protein ids to match
    :return: a dataframe containing records whose name is in the specified list
    """
    return df_utils.filter_dataframe(prot_df, col_name='protein_id', values=prot_ids)


def extract_columns(df, columns):
    """
    Create and return a dataframe containing all rows from the specified dataframe
    containing columns that match the specified column list.

    :param df: the original dataframe
    :param columns: list of columns to return
    :return: a dataframe with all rows from the original dataframe and columns
    that match the specified columns list
    """
    return df[columns]
    # common_cols = list(set(columns).intersection(set(df.columns)))
    # return df[common_cols]


def prepare_data(raw_prot_data_dir, prepared_data_dir, rna_counts_file,
                 prepared_rna_counts_file, prepared_prot_counts_file):
    """
    This method prepares the proteomic data from the raw data, by aggregating
    patients data from all individual proteomic files into a single file and
    reconciling patient ids from the protein count and RNA counts. The results
    are then saved in the specified "prepared" files. The following are the
    steps of the aggregation process:

        1) From the proteomic data file, extract the rows that contain only the
           protein count for each patient.
        2) Reconcile the patient ids in both the proteomic data and RNA data.
        3) Create an aggregate file that contains all proteomic data for all patients
           that appear in both the proteomic data and RNA data.
        4) Create an RNA count file for patient that appear in both proteomic and RNA file.
        5) Normalize the counts in both files.

    :param raw_prot_data_dir: the directory containing individual proteomic files
    :param prepared_data_dir: directory to save the prepared data
    :param rna_counts_file: the file containing the rna count
    :param prepared_rna_counts_file: the file to save the prepared and unified rna counts
    :param prepared_prot_counts_file: the  file to save the prepared and unified rna counts
    :return: a dataframe for the prepared protein count, a dataframe for the prepared rna
             counts, total number of prepared records, number of records containing any NA values
    """
    individual_prot_dir = os.path.join(prepared_data_dir, "individual_proteomic")
    utils.ensure_dir_exists(individual_prot_dir)

    # clean up the proteomic data file and save them individually
    total_cnt, na_cnt = filter_protein_data(input_dir=raw_prot_data_dir, output_dir=individual_prot_dir)

    print(f'Aggregating all proteomic data into a single file ...')
    prot_df = aggregate_proteomic_data(individual_prot_dir)

    # load the rna counts in a df
    print(f'Unifying patient records from proteomic and rna files ... ')
    rna_df = pd.read_csv(rna_counts_file, sep=FIELD_SEPARATOR)

    # make dataframes from the original ones that contains only patient data
    # that appear in both. If a patient appears in one but not the other,
    # it will be removed.
    prot_df, rna_df = unify_prot_and_rna_data(prot_df, rna_df)

    # write the new data
    print(f'Persisting proteomic data {prot_df.shape} - {[prepared_prot_counts_file]} ...')
    prot_df.to_csv(prepared_prot_counts_file, sep=COMMA, index=False, header=True)

    # write the new data
    print(f'Persisting RNA data {rna_df.shape} - {[prepared_rna_counts_file]} ... ')
    rna_df.to_csv(prepared_rna_counts_file, sep=COMMA, index=False, header=True)
    return prot_df, rna_df, total_cnt, na_cnt


def unify_prot_and_rna_data(df1, df2, normalize_col_names=True):
    """
    This function creates two new dataframes one from the protein counts and one
    from rna counts. The new dataframes contain only columns (i.e., patients) that
    exist in both dataframes. If normalize_col_names is set to ture, before choosing
    the common patient ids, the ids are normalized in the two files. Normalization
    converts the period character (.) in the ids to underscore (_) and eliminates
    the string ".T" from the end of all ids.

    :param df1: the first dataframe
    :param df2: the second dataframe
    :param normalize_col_names: indicates whether the name of the columns should
    be normalized/unified.
    :return: Copies of the original dataframes with columns that appear in both the
    original dataframes.
    """
    if normalize_col_names:
        df1.columns = utils.normalize_col_names(df1.columns.tolist())
        df2.columns = utils.normalize_col_names(df2.columns.tolist())

    # find the common columns (patient ids) from the two dataframes
    common_cols = list(set(df1.columns.tolist()).intersection(set(df2.columns.tolist())))

    # make sure column names (sample ids) are in the same order
    # so the data from every row from the two file match the sample ids
    common_cols = sorted(common_cols)

    # Extract a df with only patient ids that exist in both protein and rna files
    # add the protein id to the column list for the protein df
    df1_columns_cc = list(common_cols)
    # column zero is not patient id but the protein id
    df1_columns_cc.insert(0, df1.columns[0])
    df1 = df1.loc[:, df1_columns_cc]

    # add the rna name to the column list for the rna df
    df2_columns_cc = list(common_cols)
    df2_columns_cc.insert(0, df2.columns[0])
    df2 = df2.loc[:, df2_columns_cc]

    return df1, df2


def aggregate_proteomic_data(data_dir):
    """
    Aggregate all proteomic data into a single dataframe from all files in the specified directory.
    :param data_dir: directory containing all csv files
    :return: a datafrome containing data from all individual files
    """
    data_files = utils.get_files_in_dir(data_dir)

    df_final = None
    file_number = 0
    for data_file in data_files:
        print(f'Aggregating {os.path.basename(data_file)} ...')
        file_number += 1
        df = pd.read_csv(data_file, header="infer")

        if file_number == 1:
            df_final = df
        else:
            # we want the protein ids that appear in all files, so we use "inner"
            df_final = pd.merge(df_final, df, on='protein_id', how='inner')

    return df_final


def filter_data_file(input_file, filtered_file, separator, no_na_filtered_file=None):
    """
    The proteomic data file contain all sorts of information in them. For this project,
    we only need the count information for each protein for every patient. This method
    extract those data from the input file and stores them in the specified output file.
    The first column in the output file is the protein id and all other column contain
    patient ids as header and then the count for the protein as data on each row.

    :param input_file: path to the datafile to be filtered
    :param filtered_file: path to a file where the results should be written
    :param separator: field separator, e.g. tab or comma
    :param no_na_filtered_file: path a file where the results containing no N/A values is written
    :return: number of rows with N/A values and total rows
    """
    with open(input_file, 'r', encoding='utf-8') as source, open(filtered_file, 'w', encoding='utf-8') as target:
        if no_na_filtered_file is not None:
            no_na_target = open(no_na_filtered_file, 'w', encoding='utf-8')
        else:
            no_na_target = None

        #   The content of the file is as follows:
        #   First line contains a hash sign followed by a floating point number
        #   (looks like a version) - like #1.3
        #
        #   The second line contains information about the samples and where in the file
        #   they can be found. It contains four number, separated by a tab, as follows:
        #
        #   A       B       C       D
        #   A:  Number of RNA rows
        #   B:  Number of samples
        #   C:  The index of the last columns containing data that are not paient RNA count.
        #       Therefore, the first patients data appear at column C + 2 (starting at 1), or
        #       index C + 1 (zero based-index) - However, this is not consistent in all files.
        #       So, for the purpose of processing, we don't use this and look for a line that
        #       starts with "Sample.ID".
        #   D:  It's the number of the line where information about the patient's data appear.
        #       The First RNA symbol seem to appear 4 lines after this line. So, the first
        #       RNA symbol can be found at line D + 4.
        na_count = -1
        total_count = 0
        line_no = 0
        sample_line_found = False
        for line in source:
            line_no += 1
            # the fist line is the version line
            if line_no == 1:
                continue

            # second line has the information about the organization of the file
            if line_no == 2:
                fields = line.split(FIELD_SEPARATOR)
                first_patient_col = int(fields[2]) + 1
                first_rna_line = int(fields[3]) + 4
                continue

            # Create a header line. Information is found at the patient_id_line
            if not sample_line_found and line.startswith("Sample.ID"):
                values = utils.eliminate_fields(line, 1, first_patient_col, FIELD_SEPARATOR)
                values[0] = 'protein_id'
                header_line = separator.join(values)
                target.write(header_line)
                sample_line_found = True
                continue

            # skip lines until the first RNA line is seen
            if line_no < first_rna_line:
                continue

            # separate the fields and then pick the columns that have patient info
            values = utils.eliminate_fields(line, 1, first_patient_col, FIELD_SEPARATOR)
            new_line = separator.join(values)
            target.write(new_line)

            total_count += 1
            # add it to the na-na file only if the line contains NA for any values
            if no_na_target is not None and 'NA' in values:
                na_count += 1
                no_na_target.write(new_line)

            if 'NA' not in values:
                na_count += 1

    source.close()
    target.close()
    if no_na_target is not None:
        no_na_target.close()

    return na_count, total_count


def filter_protein_data(input_dir, output_dir):
    """
    The proteomic data files contain all sorts of information that we don't need.
    This method, process all the data in the specified "input_dir" directory whose
    name start with "hg" and end in "gct" and stores the results in the
    "output_data_dir". The prepared data contain only the protein ids (on the first
    column) and patient ids and their relevant counts for each protein id in all
    other column.

    :param input_dir: directory containing all individual proteomic data
    :param output_dir: directory where the individual processed files are to be stored
    :return: total number of records processed and number of records with NA values
    """
    print(f'Extracting patients\'s protein counts data from proteome files ')
    # make sure that the destination directory where the processed files
    # are stored exist
    utils.ensure_dir_exists(output_dir)

    # match relevant files by keyword in the name of the file
    data_type = 'proteome'
    pattern = f'^hg.*-{data_type}-.*\\.gct$'

    # get a list of all matching file name/path
    matched_file_paths = utils.get_files_in_dir(input_dir, pattern)
    na_cnt = -1
    total_cnt = 0
    for data_file in matched_file_paths:
        # filtered files get rid of all rows that are not RNA counts - data files start
        # with rows that contain information other than the counts. For this purpose,
        # we don't need them. We eliminate them before loading the data in data frames.
        f_name, _ = utils.get_file_name_and_ext(data_file)
        filtered_file = os.path.join(output_dir, f'{f_name}.csv')

        # no_na_filtered_file = os.path.join(output_dir, f'{f_name}-nona.csv')
        no_na_filtered_file = None

        na, total = filter_data_file(input_file=data_file, filtered_file=filtered_file,
                                     no_na_filtered_file=no_na_filtered_file, separator=COMMA)

        print(f'Processed \"{os.path.basename(data_file)}\" -> Total rows extracted: {total}')

        na_cnt += na
        total_cnt += total

    if na_cnt < 0:
        na_cnt = -1
    return total_cnt, na_cnt


# 1. Create a mapping for RNA Name to RNA ID, using var_map.tsv, which has the following contant:
#       RNA ID                  Description (RNA Name)      biotype
#       ENSG00000223972.5       DDX11L1                     pseudogene
#
#       Result:
#           RNA Name -> RNA ID
#           A4GALT -> ENSG00000128274.17
def create_rna_name_id_map(var_map_file):
    """
     Create a mapping for RNA Name to RNA ID, using var_map.tsv.
    :param var_map_file: the original var map containing rows in the form of, "RNA ID   Description (RNA Name)  biotype"
    :return:a map of RNA Names to RNA IDsm e.g., A4GALT -> ENSG00000128274.17
    """
    with open(var_map_file, 'r', encoding='utf-8') as input:
        rna_name_id_map = dict()
        dups = set()
        for line in input:
            tokens = line.split(FIELD_SEPARATOR)
            if tokens[1] not in dups:
                rna_name_id_map[tokens[1]] = tokens[0]
            else:
                dups.add(tokens[1])
                print(f'Found duplicate rna name: {tokens[1]}')
    return rna_name_id_map


def find_normalization_factor(rna_file, cell_rna_total_count_map):
    """
    Calculate the normalization factor based on the RNA counts from the proteomic samples
    and single cell RNAs. The calculation is done according to the following:

    For each RNA, r_i, calculate the total across all samples (RNA count file) for their log_10 values
        - use zero for the log_10 of zero counts
        - calculate the correction factor for each RNA (r_i), as cor_fact<c, r_i> = total_r_i / total_c
    Calculate the average of all correction_factors for cell c, as avg_cor_factor<c>

    :param rna_file: path to the RNA count file
    :param cell_rna_total_count_map: a map containing the total count for each single cell id
    :return: the normalization factor
    """
    return calc_total_count_normalization_factor(pd.read_csv(rna_file, header='infer'), cell_rna_total_count_map)


def calc_total_count_normalization_factor(rna_df, cell_rna_total_count_map):
    """
    Calculate the normalization factor based on the RNA counts from the proteomic samples
    and single cell RNAs. The calculation is done according to the following:

    For each RNA, r_i, calculate the total across all samples (RNA count file) for their log_10 values
        - use zero for the log_10 of zero counts
        - calculate the correction factor for each RNA (r_i), as cor_fact<c, r_i> = total_r_i / total_c
    Calculate the average of all correction_factors for cell c, as avg_cor_factor<c>

    :param rna_df: dataframe containing the RNA counts for the protemoic measures
    :param cell_rna_total_count_map: a map containing the total count for each single cell id
    :return: the normalization factor
    """

    # a map to hold normalization factor for each cell
    single_cell_normalization_factors = {}

    i = 0
    count = len(cell_rna_total_count_map)
    for cell_id in cell_rna_total_count_map.keys():
        # for testing, just calculate for a few
        i += 1

        # For each RNA, calculate the total of log10 values across all samples total_r_i,
        # log10 is make the rna values comparable with the value from the single cell file,
        # which are in log10.
        replace_inf_with = 0.0
        rna_df_log_normalized = df_utils.log_normalize_df(rna_df, replace_inf_with)

        # the first column is the RNA name, so skip it
        norm_factors = []
        for sample_id in list(rna_df_log_normalized.columns[1:]):
            total_r = sum(rna_df_log_normalized[sample_id])
            norm_factors.append(total_r / cell_rna_total_count_map[cell_id])

        # add the mean of all normalization factors for the cell to the map
        single_cell_normalization_factors[cell_id] = stats.mean(norm_factors)
        print(f'({i} / {count}) | normalization_factor for cell {cell_id} : {norm_factors[-1]}')

    # take the average of all normalization factors for all cells
    total_avg_corr_factor = stats.mean(single_cell_normalization_factors.values())
    print(f'total_avg_corr_factor: {total_avg_corr_factor}')
    return total_avg_corr_factor


def extract_single_cell_non_zero(var_map_file, cell_rna_id_counts_file, normalization_method):
    """
    normalized_humanized_mat.tsv:
      rna_id  cell_id_1   ...     cell_id_n
      RNA1    val_11      ...     val_1n
      RNA2    val_21      ...     val_2n
      ...     ...         ...     ...
      RNA_m   val_m_1     ...     val_mn

    -----------------------------------------------------------------------------
    var_map.tsv (3 fields)
      RNA Name                Description (id)        biotype
      ENSG00000223972.5       DDX11L1                 pseudogene

    -----------------------------------------------------------------------------
    var_map_full.tsv (11 fields)
    id            id.description                  geneSymbol  protein_mw      ...
    NP_958782.1   plectin isoform 1 GN=PLEC       PLEC        533778.0        ...

    1. Match the rna_id from normalized_humanized_mat to geneSymbol from var_map_full file
    2. for each cell_id from normalized_humanized_mat, find all rna_ids that are not zero
    3. write the results if a file as:
          cell_id                             non_zero_rna_values
          BPK.12x.4NQO_AAACCTGCACCCAGTG.1     rna_1=1.3,rna_2:3.0,rna_x:1.11

    ========================================================================================
    1. Create a mapping for RNA Name to RNA ID, using var_map.tsv, which has the following contant:
          RNA ID                  Description (RNA Name)      biotype
          ENSG00000223972.5       DDX11L1                     pseudogene

          Result:
              RNA Name -> RNA ID
              A4GALT -> ENSG00000128274.17

    2. From normalized_humanized_mat.tsv, for each Cell Id, extract a list of RNA Names and
       their values for all non-zero values. Replace RNA Names with RNA IDs

      rna_name        cell_id_1     cell_id_2   ... cell_id_n
      RNA_name_1      val_1_1       val_1_2         val_1_n

      Output:
          cell_1 -> { rna_id_1: n1, ..., rna_id_m: nm }

    :param var_map_file: path to the var map file
    :param cell_rna_id_counts_file: path to the file containing the single single cell rna count
    :param normalization_method: the normalization method to use
    :return: total normalization factor
    """
    if normalization_method.lower() not in ['total', 'none']:
        raise ValueError(f'Unsupported normalization {normalization_method}. Only "total" and "none" are supported.')

    # This maps the RNA names to RNA ids
    rna_name_id_map = create_rna_name_id_map(var_map_file)

    # load the single cell rna id count in a dataframe
    single_cell_df = pd.read_csv(cell_rna_id_counts_file, sep=FIELD_SEPARATOR, header='infer')

    # rna_id column doesn't have a name in the file, rename it to rna_id
    id_col_name = 'rna_name'
    single_cell_df.rename(columns={'Unnamed: 0':id_col_name}, inplace=True)

    # this is the overall map whose keys are the cell names and whose values are maps
    # that has each RNA as key and its count as value.
    cell_rna_counts_map = dict()

    # this is to save the total RNA counts for each cell. The key is the cell id
    # and the value is the sum of all RNA that have a matching RNA name in the
    # var map file and have a non-ero value. This is used for normalization later,
    # when predicting protein values for single cells.
    cell_rna_total_count_map = dict()

    # the first column is the label for the RNA name column.
    # cell names start at column 1
    cell_columns = single_cell_df.columns[1:]

    # keeps track of all RNA ids from the single cell file that do not have
    # corresponding RNA id in the var_map. The non-matching names are skipped.
    no_match_rna_names = set()

    num_cells = len(cell_columns) - 1
    # process each column, i.e. cell id
    for col_idx in range(num_cells):
        print(f'Processing column {col_idx + 1} of {num_cells} - Cell: {cell_columns[col_idx]}')
        df = single_cell_df[[id_col_name, cell_columns[col_idx]]]

        # eliminate zero values
        df = df[df.iloc[:, 1] > 0]

        # a map to hold RNA IDs and their counts for each cell
        rna_val_map = dict()

        cell_total_rna = 0

        # each row is an rna for the cell that is being processed
        for row in range(0, df.shape[0]):
            rna_name = df.iat[row, 0]
            # if there's a matching id for the rna_name, add it to the map
            if rna_name in rna_name_id_map:
                rna_val = float(df.iat[row, 1])
                cell_total_rna += rna_val
                rna_id = rna_name_id_map[rna_name]
                # rna id -> rna value
                rna_val_map[rna_id] = rna_val
            # no matching id is found for the rna_name
            else:
                no_match_rna_names.add(rna_name)

        # add the map for the cell to the overall map
        cell_rna_counts_map[cell_columns[col_idx]] = rna_val_map

        # add the total rna count for the cell to the total map
        cell_rna_total_count_map[cell_columns[col_idx]] = cell_total_rna

    print(f'Did not find RNA Id matched for: {no_match_rna_names}')
    # return
    #   - a map whose key is the single cell id and whole value is map of rna ids and their counts
    #   - a set containing rna names that do not have a matching value in the var map
    #   - a map of cell ids and total for all RNAs for that cell
    return cell_rna_counts_map, no_match_rna_names, cell_rna_total_count_map


def load_single_cell_counts(single_cell_counts_json_file):
    """
    load a map of cell_id to RNA count from the specified file that contains the
    json representation of the map. The json content must be in the following
    format.
      {
          cell_id: {
              rna_id_1: count_1
              ...
              rna_id_n: count_n
          }
      }
    :param single_cell_counts_json_file: path single cell rna count
    :return: a map constructed from the content of the file (json)
    """
    with open(single_cell_counts_json_file, 'r', encoding='UTF-8') as input_json:
        return json.load(input_json)
