# RNAToProtein
A linear regression framework for predicting proteomic values from transcriptomic data.

VERVIEW
--------
The goal of this project is to predict protein values for single cells based on their RNA counts.
The program works as follows by building a model from RNA and Protein samples measure for
patients. It extracts all samples for the specified protein from the Protein and RNA samples,
builds a linear regression model and uses it to predict the protein value for a single cell
based on its RNA values. The model simply takes the RNA values and produces a protein value.

It is assumed that the counts for single cells are log normalized but the numbers in the
RNA counts from the proteomic data sets are not. Thus, the RNA counts from the proteomic
data sets are converted to log based 10 before making and training regression models.

It is further assumed that the RNA counts are proportional to the total counts. Because
the total counts from the proteomic and single cell RNA don't match, we calculate an
adjustment factor and multiply the single cell counts with it to try to address this
issue. This factor is called the "Normalization factor" in the code and can be calculated
using the calc-normalization-factor command of the program.

HOW TO EXECUTE
--------------
The program supports five different commands. The following is how the program is executed:

    prot_predictor  <command>   <options>

Supported commands and options:

    prepare
        Reads all raw data stored in a directory whose name matches a certain format,
        eliminates all rows and columns that are not relevant and saves the results
        in individual files. The content of these files has the RNA id in the first
        columns and the counts for each sample in the following columns.

        Command:
                prepare
        Options:
            --config
                The file containing the program's global configuration values
        Example:
            python3 prot_predictor prepare \
                --config /Users/joe/config.json

    --------------------------------------------------------------------------------

    extract
        This command take a list of RNA ids, pulls rows with matching ids from all
        prepared files and aggregates them into a file specified in the output-path.

        Command:
            extract
        Options:
            --prot_ids
                list of protein ids for which RNA samples should be extract
            --output-path
                path to a directory under which the results should be saved
            --config
                The file containing the program's global configuration values

        Example:
            python3 prot_predictor extract \
                --config /Users/joe/config.json \
                --prot_ids NP_892006.3 NP_958782.1 \
                --output_dir /tmp/output_1

    --------------------------------------------------------------------------------

    model
        Make and Test linear regression models for each proteinbased on RNA values from tumor counts.
        The result is a CSV file with the protein id, model's properties, prediction error and r2 score.

        Command:
            model
        Options:
            --prot_ids
                Optional list of protein ids. If omitted, the program runs for all proteins in the original file
            --config
                The file containing the program's global configuration values
        Example:
            python3 prot_predictor model \
                --config ~/config.json \
                --prot_ids NP_892006.3  NP_958782.1 NP_002464.1

    --------------------------------------------------------------------------------

    calc-normalization-factor
        Calculate the normalization factor and write the value in a file to be used later
        for prediction. This value ise used to normalize the rna values for single cells
        to be comparable to RNA values from the proteomic data.

        Command:
            calc-normalization-factor
        Options:
            --config
                The file containing the program's global configuration values
        Example:
            python3 prot_predictor calc-normalization-factor \
                --config /Users/joe/config.json

    --------------------------------------------------------------------------------

    prep-single-cell
        Prepare Single Cell RNA values to be used for protein value prediction. This is a prerequisite step
        to predicting any protein value based on single cell RNA values.

        Command:
            prep-single-cell
        Options:
            --config
                The file containing the program's global configuration values
        Example:
            python3 prot_predictor prep-single-cell \
                --config /Users/joe/config.json \

    --------------------------------------------------------------------------------

    predict-single-cell
        Predict the protein values for single cells based on their RNA measurements

        Command:
            predict-single-cell
        Options:
            --prot_ids
                list of proteins to predict for a cell
            --cell_ids
                the id of a cell to predict for
            --output_dir
                optional name for directory to store the results. The directory is created under the "generated"
                directory as specified in the config file. If this option is omitted, results are printed on screen
            --config
                The file containing the program's global configuration values
        Example:
            python prot_predictor predict-single-cell \
                --config /Users/joe/config.json \
                --prot_ids  NP_000436.2 NP_112598.3 NP_001333374.1 \
                --cell_ids BPK.12x.4NQO_AAACCTGCACCCAGTG.1 \
                --output_dir pred_3



REMARKS:
    A config file must be supplied to the program for all commands. The config file contains information
    that is necessary for the program to operate. You can control some behavior of the program by setting
    or changing the configuration values in the configuration file.

    The first step before running any other command is to run the "prepare". The prepare command
    processes the raw tumor data and produces aggregate data files that will be used by all other
    command.

    Before attemtping to predict a single cell protein value, you must have run "prepare",
    "prepare-single-cell", and "calc-normalization-factor". Once you run those commands in the
    sequence, you can run the prediction command as many times without having to run those other
    commands, unless the underlying data changes.


Normalization Factor calculation method (according to professor Babur's instructions):
    - We will rescale the single-cell RNA values by assuming that the total RNA is not changing from sample to sample.

    - The single-cell data has a lot of missing values, so we cannot expect the non-missing values will add
      up to the same total. They will add up to a much lower value, even after a proper normalization, because
      many values are missing.

    - Let's focus on the non-missing values in a cell (say cell_1, the first column in the single cell data).
      Let's assume the total of these RNA will not change from cell to cell or from sample to sample.
      Now we can look at the corresponding values in one of the training samples (say sample_1), calculate
      the total, and come up with a correction factor for cell_1 to make the two totals equal.

      When we do this with sample_2 instead, we will probably find another correction factor for cell_1.
      Let's find all  the correction factors (corresponding to each sample), and use their average for cell_1.
      Let's call this cell_1_cor_fac.

    - Do the same thing for cell_2 and find cell_2_cor_fac, and continue doing for all cells.

    - Having a different correction factor for each cell does not make sense because cells are already normalized
      among themselves. Instead, let's average all these correction factors cell_1_cor_fac, cell_2_cor_fac,
      cell_3_cor_fac, and such, and generate that one correction factor (say cor_fac) to use for all the cells.

    - To normalize the single cell data (and make it comparable with our training data), we will simply multiply
      each value in the matrix with the cor_fac.
