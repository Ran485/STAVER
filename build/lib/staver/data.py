import os
import pandas as pd
import numpy as np

# import dask.dataframe as dd
from pyteomics import mzid, pepxml, protxml, mztab
from utils import timer, status_info, reduce_mem_usage
from joblib import Memory, Parallel, delayed
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from rich.progress import Progress, track, BarColumn
from rich.console import Console


# Get the number of CPU threads
NUM_WORKERS = cpu_count() - 2
# print(f'Calling {NUM_WORKERS} CPU threads for parallel processing.')


def find_all_files(input_folder, extension=".tsv"):
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(extension):
                fullname = os.path.join(root, file)
                # print(fullname)
                yield fullname


def get_all_files(input_folder, extension=".tsv"):
    all_filenames = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(extension):
                fullname = os.path.join(root, file)
                all_filenames.append(fullname)
    return all_filenames


@timer
@status_info()
def pool_load_file(inpath=None, extension=".tsv"):
    all_filenames = find_all_files(inpath, extension)
    combined_csv = pd.concat(map(read_file, all_filenames))
    return combined_csv


def read_proteomics_file(file_path, usecols):
    """
    Reads various proteomics output file formats and returns the content.

    This function reads data from different proteomics file formats and returns
    it as a pandas DataFrame. It supports common file formats like CSV, TSV, Excel,
    and specialized formats like mzid, pepxml, protxml, and mztab.

    Parameters:
        file_path (str): Path to the proteomics file.

    Returns:
        pd.DataFrame: DataFrame containing the data from the file.

    Raises:
        ValueError: If the file does not exist or if the file extension is unsupported.
    """
    if not os.path.exists(file_path):
        raise ValueError(f"File {file_path} does not exist!")

    ext = os.path.splitext(file_path)[1].lower()

    # Reading common file formats
    if ext in [".txt", ".tsv", ".csv"]:
        delimiter = "\t" if ext == ".tsv" else ","
        return pd.read_csv(file_path, usecols=usecols, delimiter=delimiter)

    elif ext in [".xls", ".xlsx", ".xlsm"]:
        return pd.read_excel(file_path, usecols=usecols)

    # Reading specialized proteomics file formats
    try:
        if ext == ".mzid":
            return pd.DataFrame(mzid.read(file_path))

        elif ext == ".pepxml":
            return pd.DataFrame(pepxml.read(file_path))

        elif ext == ".protxml":
            return pd.DataFrame(protxml.read(file_path))

        elif ext == ".mztab":
            with mztab.MzTab(file_path) as reader:
                psm_data = reader.spectrum_match_table
                return pd.DataFrame(psm_data)

    except Exception as e:
        raise ValueError(f"Error reading {ext} file at {file_path}: {str(e)}") from e

    raise ValueError(f"Unsupported file extension: {ext}")


def read_file(file_path, filter_protein_qvalue=False, verbose=False):
    """
    Reads a TSV file and returns a filtered pandas DataFrame.

    This function reads a TSV file into a DataFrame, optionally filters rows based on Q.Value
    and Protein.Q.Value thresholds, and adds a column with the file name. It also optimizes
    memory usage of the DataFrame.

    Args:
        file_path (str): The path to the TSV file.
        filter_protein_qvalue (bool, optional): If True, filters the DataFrame based on both
                                                Q.Value and Protein.Q.Value being less than 0.05.
                                                If False, filters based only on Q.Value. Defaults to False.
        verbose (bool, optional): If True, prints error messages. Defaults to False.

    Returns:
        pandas.DataFrame or None: The DataFrame containing the TSV file data, or None if an error occurs.

    Raises:
        Exception: If an error occurs during file reading or processing.
    """
    try:
        # Define the columns to be read
        columns_to_read = [
            "File.Name",
            "Genes",
            "Stripped.Sequence",
            "Precursor.Id",
            "Modified.Sequence",
            "Precursor.Normalised",
            "Q.Value",
            "Protein.Q.Value",
        ]

        # Read the TSV file
        df = read_proteomics_file(file_path, usecols=columns_to_read)

        # Filter the DataFrame based on Q values
        if filter_protein_qvalue:
            df = df[(df["Q.Value"] < 0.05) & (df["Protein.Q.Value"] < 0.05)]
        else:
            df = df[df["Q.Value"] < 0.05]

        # Optimize memory usage
        df = reduce_mem_usage(df, verbose=False)

        # Add the file name column
        df["file_name"] = file_path.split("/")[-1]

        return df

    except Exception as e:
        if verbose:
            print(f"Error while processing `{file_path}`: {e}")
        return None  # Return None explicitly if there's an error


# @timer
def ThreadPool_load_file(inpath=None, extension=".tsv"):
    all_filenames = get_all_files(inpath, extension)
    print(f"\nCalling {NUM_WORKERS} CPU threads for parallel processing.")
    with ThreadPoolExecutor(NUM_WORKERS) as pool, Progress(
        BarColumn(bar_width=None, complete_style="green", finished_style="green"),
    ) as progress:
        task_list = [
            progress.add_task("Processing", filename=file.split("/")[-1])
            for file in all_filenames
        ]
        results = list(pool.map(read_file, all_filenames))
        for task in task_list:
            progress.update(task, completed=1)
    return pd.concat(results)


def joblib_load_file_track(inpath=None, extension=".tsv"):
    all_filenames = find_all_files(inpath, extension)
    print(f"\nCalling {NUM_WORKERS} CPU threads for parallel processing.")
    with Parallel(n_jobs=NUM_WORKERS) as parallel:  # 使用多进程模式
        results = parallel(
            delayed(read_file)(file)
            for file in track(all_filenames, description="Parallel Processing")
        )
    return pd.concat(results)  # 合并所有的DataFrame


# @timer
# @status_info()
def joblib_load_file(inpath=None, extension=".tsv", num_workers=NUM_WORKERS):
    """
    Loads multiple files in parallel using joblib and returns a concatenated DataFrame.

    Args:
        inpath (str): The path to the directory containing the files. Defaults to None.
        extension (str): The file extension to filter the files. Defaults to ".tsv".
        num_workers (int): The number of CPU threads for parallel processing. Defaults to NUM_WORKERS.

    Returns:
        pandas.DataFrame: The concatenated DataFrame.

    Examples:
        >>> joblib_load_file(inpath="/path/to/files", extension=".csv", num_workers=4)
            Column1  Column2
        0        1        2
        1        3        4
    """
    all_filenames = get_all_files(inpath, extension)
    console = Console()
    console.print(
        f"\nCalling {num_workers} CPU threads for parallel processing...",
        style="bold yellow",
    )

    results = []
    with Progress() as progress:
        task_id = progress.add_task(
            "[cyan]Parallel loading files...", total=len(all_filenames)
        )

        with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
            futures = {executor.submit(read_file, file): file for file in all_filenames}

            for future in as_completed(futures):
                results.append(future.result())
                progress.update(task_id, advance=1)

    return pd.concat(results)


if __name__ == "__main__":
    DIA_PATH = "/Users/ranpeng/Desktop/DIA-QC/data/likai-diann-raw-20/"
    data = joblib_load_file(inpath=DIA_PATH)
    print(data.shape)
    print(data.head())
