import os
from os import walk
import re
import zipfile
import datetime


#
# This file contains a few utility functions that are used for file and directory
# related operations.
#
# @author: Cyrus
#

def get_files_in_dir(root_dir, pattern=None):
    """
    Get a list of all files that are in any of the sub-folders of the specified root directory.
    If the "pattern" parameter is not None, only return the files whose name matched the pattern.

    :param root_dir: the root directory
    :param pattern: the pattern to which the name of the files should be matched
    :return: a list containing all files in any of the sub-folders of the specified directory,
    whose name match the specified pattern.
    """
    if pattern is not None:
        type_pattern = re.compile(pattern)
    else:
        type_pattern = None

    file_paths = list()
    w = walk(root_dir)
    for(dir_path, dir_names, file_names) in w:
        for file_name in file_names:
            if type_pattern is None or not re.search(type_pattern, file_name):
                file_paths.append(os.path.join(dir_path, file_name))

    return file_paths


def get_file_base_name(file_path):
    """
    Get the base name of the specified file path

    :param file_path: file path
    :return: the base name from the file path
    """
    ridx = file_path.rindex(".")
    ext = file_path[ridx + 1:]
    base = file_path[0: ridx]
    return base, ext


def get_file_name_and_ext(file_path):
    """
    Get the file name and its extension from the specified path.

    :param file_path: full path to the file
    :return: file name and its extension
    """
    filename = os.path.split(file_path)[-1]
    return os.path.splitext(filename)


def ensure_dir_exists(directory_path):
    """
    Make sure that specified directory exits. If the specified directory exists, this
    function doesn't do anything. If it doesn't exist, it created is.

    :param directory_path: the path to the directory to check
    :return: None
    """
    if not os.path.isdir(directory_path):
        os.makedirs(directory_path)


def extract_field_values(line, from_idx, to_idx, separator):
    """
    This method splits the specified line into a list of values using specified separator
    and returns a sub-list starting at from_idx (inclusive) to to_idx (exclusive).

    :param line: string to parse and extract a sub-list from
    :param from_idx: the starting index (inclusive)
    :param to_idx: the ending index (exclusive)
    :param separator: field separator
    :return: list containing specified elements from the line if indices are valid/exist. An empty list otherwise
    """
    fields = line.split(separator)
    if 0 <= from_idx <= to_idx and len(fields) < to_idx:
        return fields[from_idx : to_idx]

    return []


def eliminate_fields(line, from_idx, to_idx, separator):
    """
    Parse the line into tokens using the specified separator and then eliminate the values
    starting at index "from_idx" up to but not including "to_idx". So, a call to this method
    with the following parameters:
      line = 'a,b,c,d,e'
      separator = ','
      from_idx = 1
      to_idx = 3
    will return [a, d, e]

    :param line: a string containing values
    :param from_idx: the starting index (inclusive)
    :param to_idx: the ending index (exclusive)
    :param separator: field separator
    :param separator:
    :return: a list with target values from the original line
    """
    fields = line.split(separator)
    if 0 <= from_idx <= to_idx and len(fields) < to_idx:
        return fields[from_idx : to_idx]
    # if there's anything that should be saved from the beginning of the list
    if from_idx > 0:
        vals = fields[0: from_idx]
    else:
        vals = []
    # if there's anything that should be saved from the end of the list
    if to_idx < len(fields):
        return vals + fields[to_idx:]
    return vals


#
#
def normalize_col_names(columns):
    """
    Modify column names as follows:
      - Replace the dot character
      - remove string ".T" from the end of column names

    :param columns: column names
    :return: modified column names
    """
    new_cols = []
    for s in columns:
        new_cols.append(re.sub('\\.', '_', re.sub('\\.T$', '', s)))
    return new_cols


def unzip(zip_file, output_dir):
    """
    Unzip the specified file into the specified folder.

    :param zip_file: tha path to the file to be unzipped
    :param output_dir: the path to the folder to write the results
    :return: None
    """
    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
        zip_ref.extractall(output_dir)


def convert_time_in_hours_mins(secs):
    """
    Converts the specified seconds to number of hours and number of minuts.

    :param secs: the number of seconds to convert
    :return: a human readable string representing the number of seconds in term of hours and minutes
    """
    # if secs == 0:
    #     str = '0 second'

    hours = secs // 3600
    minutes = (secs % 3600) // 60
    secs = secs % 60
    str = ''
    if hours > 0:
        str = f'{hours} hour(s)'
    if minutes > 0:
        sep = ',' if len(str) > 0 else ''
        str += f'{sep} {minutes} minutes(s)'
    # if secs > 0:
    sep = ',' if len(str) > 0 else ''
    str += f'{sep} {secs} second(s)'
    return str


def get_date_time_str():
    """
    Get the current time in a string formatted as year_month_day_hour_minutes_seconds

    :return: the current time in a string formatted as year_month_day_hour_minutes_seconds
    """
    now = datetime.datetime.now()
    return f'{now.year}_{now.month}_{now.day}_{now.hour}_{now.minute}_{now.second}'


def read_line_from_file(file_path, line_number):
    """
    Read and return the line at the specified line from the specified file.

    :param file_path: the path to the file
    :param line_number: the number of the line to read inside the file
    :return: the line located at the specified number/location in the file
    """
    with open(file_path, 'r', encoding='UTF-8') as input_file:
        offset = -1
        for i in range(line_number):
            if i == line_number - 1:
                return input_file.readline().strip()
            input_file.readline()
            if offset == input_file.tell():
                return None
            offset = input_file.tell()
