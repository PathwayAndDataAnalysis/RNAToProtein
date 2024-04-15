import numpy as np
import pandas as pd

FIELD_SEPARATOR = '\t'
COMMA = ','


#
# This file contains some utility functions related to Panda Data Frames.
#
# @author: Cyrus
#
def read_df(data_file, sep=COMMA):
    """
    create and return dataframe form the specified file using the specified separator
    :param data_file: path to the file containing the data
    :param sep: field separator, defaults to comma
    :return: a dataframe
    """
    return pd.read_csv(data_file, sep=sep)


def filter_dataframe(df, col_name, values=None):
    """
    Extract from the specified dataframe records whose name match any of the
    values specified in the values list. If the values parameter is None, the
    original dataframe is returned.

    :param df: dataframe containing the records
    :param col_name: the name of the column to match
    :param values: list of values to match for the specified column
    :return: a dataframe containing records whose value is in the specified list
    """
    if values is not None:
        return df[df[col_name].isin(values)]
    return df


def log_normalize_df(df, skip_cols_count, replace_inf_with=None):
    """
    Normalize the values of the columns in the df - log_10(val)

    :param df: the dataframe containing the data to be normalized
    :param skip_cols_count: number of columns to skip starting at 0
    :param replace_inf_with: the value to replace -inf values after normalization. np.inf is the result
    of taking the log_10 of zero.
    :return: A copy of the original dataframe whose values are the log based 10 of the values in the
    original dataframe and -inf values replaced by 0.
    """
    df2 = df.copy()
    for col in df.columns[skip_cols_count:]:
        df2[col] = df[col].apply(np.log10)
    if replace_inf_with is not None:
        return df2.replace(-np.inf, replace_inf_with)
    return df2
