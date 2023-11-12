#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File     : utils.py
@Time     : 2023/04/21 13:26:04
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
"""
# here put the import lib

import os
import sys
import time
import pandas as pd
import numpy as np
import config_file as cfg_file

from datetime import datetime
from functools import wraps
from rich.progress import track
from rich.console import Console
from joblib import Parallel, delayed
from rich.progress import track
from multiprocessing import cpu_count

# from plots import *


# Get the number of CPU threads
def get_cpu_count():
    """Return the number of CPU threads."""
    NUM_WORKERS = cpu_count() - 2
    print(f"Calling {NUM_WORKERS} CPU threads for parallel processing.")
    return NUM_WORKERS


def change_root_dir(path=None):
    """
    change the root directory: set as the 'PhosphoPreprocessy/'
    """
    print("Current Working Directory ", os.getcwd())
    # path = '/Desktop/CCA/Metadata/phosphosite-mapping（2021-11-11）/PhosphoPreprocessy'
    try:
        os.chdir(path)
        print("Current working directory: {0}".format(os.getcwd()))
    except FileNotFoundError:
        print("Directory: {0} does not exist".format(path))
    except NotADirectoryError:
        print("{0} is not a directory".format(path))
    except PermissionError:
        print("You do not have permissions to change to {0}".format(path))


# Create the output directory
def create_output_dir(filename=None, creat_time_subdir=False):
    """Create the output directory.

    Args:
        filename (str): A given filename.
        creat_time_subdir (bool, optional):
            creat 2021-11-12 subdirectory,defaults to True.
    Returns:
        output_dir (str): The output directory.
    """
    root_dir = "./"
    out_path = root_dir + filename
    if not os.path.isdir(out_path):  # in case root_dir doesn't exist
        os.makedirs(out_path)
        print(f"Successfully created output subdir: {out_path}")
        if creat_time_subdir:
            date_string = datetime.now().strftime("%Y_%m_%d")
            out_path = os.path.join(out_path, date_string)
            if not os.path.isdir(out_path):
                os.mkdir(out_path)
                print(f"Successfully created output subdir: {out_path}")
    else:
        print(f"The current path: {out_path} already exist")
    return out_path + "/"


# create a timer decorator
def timer(func):
    """Return the function start and finished cost time"""

    @wraps(func)
    def wrapper(*args, **kwargs):
        title_start = f" {func.__name__!r} function begins running... "
        print(f"\n{title_start.center(84, '=')}\n")
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        title_end = f" {func.__name__!r} function was succesfully done in {round((end_time-start_time), 2)}s "
        print(f"\n{title_end.center(84, '=')}\n")
        return result

    return wrapper


def status_info():
    """Return a status message for processing the files."""

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            console = Console()
            with console.status(
                f"[bold green]Working on running {func.__name__}...",
                spinner="aesthetic",
            ) as status:
                console.log(f"[bold][blue]Processing {func.__name__} started...")
                result = func(*args, **kwargs)
                console.log(f"[bold][red]Finished!")
                return result

        return wrapper

    return decorator


def memory():
    """Check the memory usage."""
    import psutil

    mem = psutil.virtual_memory()
    zj = float(mem.total) / 1024 / 1024 / 1024
    ysy = float(mem.used) / 1024 / 1024 / 1024
    kx = float(mem.free) / 1024 / 1024 / 1024
    print("Total system memory:%d.3GB" % zj)
    print("The system has used memory:%d.3GB" % ysy)
    print("System free memory:%d.3GB" % kx)


def reduce_mem_usage(df, verbose=True):
    # List of numeric data types
    numerics = ["int8", "int16", "int32", "int64", "float16", "float32", "float64"]
    # Starting memory usage
    start_mem = df.memory_usage().sum() / (1024**2)

    # Iterate over columns in the DataFrame
    for col in df.columns:
        # Get the dtype of the column
        col_type = df[col].dtype

        # Check if the column is numeric
        if col_type in numerics:
            # Calculate the min and max values of the column
            c_min = df[col].min()
            c_max = df[col].max()

            # Determine the optimal data type for the column
            if col_type.kind == "i":
                optimal_type = next(
                    (
                        int_type
                        for int_type in (np.int8, np.int16, np.int32, np.int64)
                        if (
                            c_min > np.iinfo(int_type).min
                            and c_max < np.iinfo(int_type).max
                        )
                    ),
                    np.int64,
                )
            else:
                optimal_type = next(
                    (
                        float_type
                        for float_type in (np.float16, np.float32, np.float64)
                        if (
                            c_min > np.finfo(float_type).min
                            and c_max < np.finfo(float_type).max
                        )
                    ),
                    np.float64,
                )
            # Convert the data type of the column to the optimal data type
            df[col] = df[col].astype(optimal_type)
    # Ending memory usage
    end_mem = df.memory_usage().sum() / (1024**2)
    # Print memory usage information
    memory_usage_info(verbose, start_mem, end_mem)
    # Return the DataFrame
    return df


def memory_usage_info(verbose, start_mem, end_mem):
    if verbose:
        print(
            "Mem. usage decreased to {:5.2f} Mb ({:.1f}% reduction)".format(
                end_mem, 100 * (start_mem - end_mem) / start_mem
            )
        )


def index_transform(df, convert_reverse=False):
    if convert_reverse:
        df["index"] = df.index
        df = df.loc[:, ~df.columns.duplicated()]
        df[["Genes", "Stripped.Sequence", "Modified.Sequence", "Precursor.Id"]] = (
            df["index"].astype(str).str.split("_", expand=True)
        )
        df.set_index(
            ["Genes", "Stripped.Sequence", "Modified.Sequence", "Precursor.Id"],
            inplace=True,
        )
    else:
        df["index"] = (
            df["Genes"]
            + "_"
            + df["Stripped.Sequence"]
            + "_"
            + df["Modified.Sequence"]
            + "_"
            + df["Precursor.Id"]
        )
        df.set_index("index", inplace=True)
    return df


def applyParallel(dfGrouped, func):
    """Apply a function to each group in a pandas dataframe in parallel.

    Args:
        dfGrouped (DataFrameGroupBy): The pandas dataframe grouped by some column.
        func (function): The function to apply to each group.

    Returns:
        DataFrame: The result of the function applied to each group.
    """
    # print(f'Calling {NUM_WORKERS} CPU threads for parallel processing.')

    # Run the function in parallel with the specified number of workers
    res = Parallel(n_jobs=NUM_WORKERS)(
        delayed(func)(group)
        for name, group in track(dfGrouped, description="Parallel processing...")
    )
    return pd.concat(res)


def save_log(outpath="./"):
    """save log to file

    Args:
        path (str, optional): The log file outpath. Defaults to './'.
    """

    class Logger(object):
        def __init__(self, filename="Default.log", outpath="./"):
            self.terminal = sys.stdout
            self.log = open(
                os.path.join(outpath, filename),
                "a",
                encoding="utf8",
            )

        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)

        def flush(self):
            pass

    fileName = datetime.now().strftime("day_" + "%Y_%m_%d")
    sys.stdout = Logger(fileName + ".log", outpath=outpath)

    # 这里输出之后的所有的输出的print 内容即将写入日志
    print(fileName.center(80, "="))


def data_transform(df, convert_reverse=False) -> pd.DataFrame:
    """
    Perform index transformation on the input DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame to be transformed.
    convert_reverse : bool, optional
        Whether to perform reverse transformation, by default False.

    Returns
    -------
    pd.DataFrame
        The transformed DataFrame.
    """
    if convert_reverse:
        df["index"] = df.index
        df = df.loc[:, ~df.columns.duplicated()]
        df[
            [
                "Genes",
                "Stripped.Sequence",
                "Modified.Sequence",
                "Precursor.Id",
                "file_name",
            ]
        ] = (
            df["index"].astype(str).str.split("_", expand=True)
        )
        df.set_index(
            [
                "Genes",
                "Stripped.Sequence",
                "Modified.Sequence",
                "Precursor.Id",
                "file_name",
            ],
            inplace=True,
        )
    else:
        df["index"] = (
            df["Genes"]
            + "_"
            + df["Stripped.Sequence"]
            + "_"
            + df["Modified.Sequence"]
            + "_"
            + df["Precursor.Id"]
            + "_"
            + df["file_name"]
        )
        df.set_index("index", inplace=True)
    return df


def median_normalize(df):
    """
    Perform median normalization on a pandas DataFrame.

    Parameters
    ----------
    df : pandas DataFrame
        pandas DataFrame to be normalized. Each column is treated as a sample, and each row is a quantity.

    Returns
    -------
        normalized pandas DataFrame
    """
    print("Performing median normalization ...\n")

    # Calculate the global median across all columns/samples
    global_median = df.median().median()

    # Subtract each column's median from all values in the column, and then add the global median.
    df_normalized = df.subtract(df.median()).add(global_median)

    return df_normalized


def split_n(s, delimiter, n):
    parts = s.split(delimiter)
    return parts[:n] + [delimiter.join(parts[n:])]


def quantile_normalize(data, method="mean"):
    """quantilen_normalize for the proteomics data
    Args:
    ------------
        df: -> Dataframe
        method: -> Default to `median`, optioanl is `mean`.

    Return:
    ------------
        Dataframe

    """
    ## Mean Quantile Normalize
    if method == "mean":
        df = data.copy()
        ranks = df.rank(method="first").stack()
        rank_mean = df.stack().groupby(ranks).mean()
        # Add interpolated values in between ranks
        finer_ranks = (rank_mean.index + 0.5).to_list() + rank_mean.index.to_list()
        rank_mean = rank_mean.reindex(finer_ranks).sort_index().interpolate()
        return df.rank(method="average").stack().map(rank_mean).unstack()
    ## Median Quantile Normalize
    elif method == "median":
        df = data.copy()
        # compute rank
        dic = {}
        for col in df:
            dic[col] = df[col].sort_values(na_position="first").values
        sorted_df = pd.DataFrame(dic)
        # rank = sorted_df.mean(axis = 1).tolist()
        rank = sorted_df.median(axis=1).tolist()
        # sort
        for col in df:
            # compute percentile rank [0,1] for each score in column
            t = df[col].rank(pct=True, method="max").values
            # replace percentile values in column with quantile normalized score
            # retrieve q_norm score using calling rank with percentile value
            df[col] = [
                np.nanpercentile(rank, i * 100) if ~np.isnan(i) else np.nan for i in t
            ]
        return df


def coefficient_variation(data=None, log_transformed=False, dropna_thresh=0.5):
    """caculate coefficient of variation for each plasma protein
    Args:
    -----------
        log_transformed -> bool: log10 transformation; default to False.
        dropna_thresh -> float: drop  plasma protein in row, default to 0.5.

    Return:
    -----------
        Dataframe

    """
    if log_transformed:
        data = np.log10(data)
    data = data.dropna(axis=0, thresh=len(data.columns) * dropna_thresh)
    data["Coefficient of Variation [%]"] = data.std(axis=1) / data.mean(axis=1)
    data["Abundance"] = data.iloc[:, :-1].mean(axis=1)
    data["Coefficient of Variation [%]"] = data["Coefficient of Variation [%]"] * 100
    data = data.sort_values("Abundance", ascending=False)
    data["Abundance_Rank"] = [i for i in range(1, len(data) + 1)]
    return data


def preprocess(data, filter_protein_qvalue=False, outpath=None):
    if filter_protein_qvalue:
        data1 = data[data["Protein.Q.Value"] < 0.05]
    data_long = data_transform(data1)
    data_long.reset_index(inplace=True)
    data_long.rename(columns={"index": "peptide"}, inplace=True)
    # Use the pivot method to convert it to wide format
    df_wide = data_long.pivot_table(
        index="peptide", columns="File.Name", values="Precursor.Normalised"
    )
    df_wide.replace(0, np.nan, inplace=True)
    df_wide = np.log2(df_wide + 1)
    df_wide_normalize = median_normalize(df_wide)
    # Melt wide format data into long format data
    df_long = df_wide_normalize.reset_index().melt(
        id_vars="peptide", var_name="File.Name", value_name="Precursor.Normalised"
    )
    df_long.dropna(inplace=True)
    df_long.rename(columns={"peptide": "index"}, inplace=True)
    df_long["str"] = df_long["index"].apply(lambda x: split_n(x, "_", 4))
    df_long[
        ["Genes", "Stripped.Sequence", "Modified.Sequence", "Precursor.Id", "file_name"]
    ] = pd.DataFrame(df_long.str.tolist(), index=df_long.index)
    df_long.drop(columns=["str", "index"], inplace=True)

    # plot the results
    print(f"Plot the density plot for raw data...\n")
    # Box and density plot before normalization
    density_plot(
        df_wide,
        title="Density curve and Rug Plot of raw data",
        x_title="LFQ Abundance [Log2]",
        y_title="Probability density",
        outpath=outpath,
    )
    print(f"Plot the boxplot for raw data...\n")
    boxplot(
        df_wide,
        log_transform=False,
        title="Box Plot of raw data",
        x_title="Samples",
        y_title="LFQ Abundance [Log2]",
        outpath=outpath,
    )

    # Box and density plot after normalization
    print(f"Plot the density plot after median normalize..\n")
    density_plot(
        df_wide_normalize,
        title="Density curve and Rug Plot after normalization",
        x_title="LFQ Abundance [Log2]",
        y_title="Probability density",
        outpath=outpath,
    )
    print(f"Plot the density plot after median normalize...\n")
    boxplot(
        df_wide_normalize,
        log_transform=False,
        title="Box Plot after median normalization",
        x_title="Samples",
        y_title="LFQ Abundance [Log2]",
        outpath=outpath,
    )

    return df_long


def measure_performance(func, sizes, data_generator):
    """
    Measures the execution time of a given function over varying sample sizes.

    This function measures the time it takes for a given function to execute with different sizes of input data,
    generated by a specified data generating function.

    Args:
        func (callable): The function whose performance is to be measured. It should accept a single argument, the data.
        sizes (list of int): A list of sample sizes for which the function's performance will be measured.
        data_generator (callable): The data generating function. It should accept the sample size as an argument and return the data.

    Returns:
        list of float: A list containing the execution times for each sample size.

    Usage:
        >>> def dummy_function(data):
        ...     time.sleep(1)  # Simulate some computation
        >>> def generate_dummy_data(size):
        ...     return [0] * size  # Generate dummy data of a given size
        >>> sample_sizes = [10, 50, 100]
        >>> times = measure_performance(dummy_function, sample_sizes, generate_dummy_data)
        >>> plt.plot(sample_sizes, times, marker='o')
        >>> plt.xlabel('Sample Size')
        >>> plt.ylabel('Time (seconds)')
        >>> plt.title('Function Performance with Varying Sample Sizes')
        >>> plt.grid(True)
        >>> plt.show()
    """
    times = []
    for size in sizes:
        data = data_generator(size)
        start_time = time.time()
        func(data)
        end_time = time.time()
        times.append(end_time - start_time)
    return times


def generate_test_peptide_data(n_peptides, n_samples, missing_data_rate=0.1):
    """
    Generates proteomics data with missing values.

    This function creates a DataFrame of proteomics data (e.g., peptide or protein abundances)
    for a given number of peptides and samples. It introduces missing values to simulate a
    typical proteomics dataset.

    Args:
        n_peptides (int): The number of peptides (or proteins), corresponding to the number of rows.
        n_samples (int): The number of samples, corresponding to the number of columns.
        missing_data_rate (float, optional): The proportion of missing values in each sample. Defaults to 0.1.

    Returns:
        pd.DataFrame: A DataFrame containing the generated proteomics data. Each row represents a peptide (or protein),
                      and each column represents a sample. Missing data are represented by NaN values.

    Example:
        >>> peptide_data = generate_peptide_data(100, 10, 0.1)
        >>> print(peptide_data.head())
    """
    np.random.seed(0)
    data = np.random.lognormal(mean=2, sigma=0.5, size=(n_peptides, n_samples))

    # Introduce missing values, ensuring approximately 'missing_data_rate' missingness per sample
    for col in range(n_samples):
        missing_indices = np.random.choice(
            n_peptides, size=int(n_peptides * missing_data_rate), replace=False
        )
        data[missing_indices, col] = np.nan

    return pd.DataFrame(data, columns=[f"Sample_{i+1}" for i in range(n_samples)])


def create_synthetic_test_data(
    num_proteins=100, num_samples=10, max_peptides_per_protein=10
):
    """
    Creates a synthetic DataFrame mimicking DIA proteomics data.

    This function generates synthetic test data for proteomics analysis. It simulates
    multiple proteins, each with a random number of peptides, and assigns random
    intensity values to these peptides across different samples. Missing values are
    introduced to simulate a realistic DIA dataset.

    Args:
        num_proteins (int): The number of proteins to simulate. Defaults to 100.
        num_samples (int): The number of samples for each peptide. Defaults to 10.
        max_peptides_per_protein (int): The maximum number of peptides per protein. Defaults to 10.

    Returns:
        pd.DataFrame: A DataFrame representing synthetic proteomics data. Rows correspond to peptides,
                      columns to samples, and values to peptide intensities.

    Example:
        >>> synthetic_data = create_synthetic_test_data(num_proteins=100, num_samples=10, max_peptides_per_protein=10)
        >>> print(synthetic_data.head())
    """
    data = []
    peptide_info = []  # Stores information about peptides and proteins

    for protein_id in range(num_proteins):
        num_peptides = np.random.randint(1, max_peptides_per_protein + 1)
        for peptide_id in range(num_peptides):
            mean_intensity = np.random.uniform(1, 3)  # Random average intensity
            sigma_intensity = np.random.uniform(0.3, 1.0)  # Random standard deviation
            peptide_intensity = np.random.lognormal(
                mean=mean_intensity, sigma=sigma_intensity, size=num_samples
            )

            # Introduce missing values more likely in low-intensity peptides
            missing_probability = np.clip(1 / peptide_intensity, 0, 1)
            missing_indices = np.random.uniform(size=num_samples) < missing_probability
            peptide_intensity[missing_indices] = np.nan

            data.append(peptide_intensity)
            peptide_info.append(f"Protein_{protein_id}_Peptide_{peptide_id}")

    # Introduce sample-specific perturbations
    sample_specific_effects = np.random.normal(loc=0, scale=0.2, size=num_samples)
    data = np.array(data) + sample_specific_effects

    # Create DataFrame
    df = pd.DataFrame(data.transpose(), columns=peptide_info)
    df.index.name = "Sample"
    return df.T
