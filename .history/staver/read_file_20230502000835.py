he


import os
import pandas as pd
import numpy as np

# import dask.dataframe as dd
from utils import timer, status_info, create_output_dir
from joblib import Memory, Parallel, delayed
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from utils import reduce_mem_usage
from rich.progress import Progress, track
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


def read_file(file, Protein_Qvalue=False):
    """
    Reads a TSV file and returns a pandas DataFrame.

    Args:
        file (str): The path to the TSV file.
        Protein_Qvalue (bool, optional): Whether to filter the DataFrame by Protein.Q.Value. Defaults to False.

    Returns:
        pandas.DataFrame: The DataFrame containing the TSV file data.
    """
    try:
        df = pd.read_table(
            file,
            usecols=[
                "File.Name",
                "Genes",
                "Stripped.Sequence",
                "Precursor.Id",
                "Modified.Sequence",
                "Precursor.Normalised",
                "Q.Value",
                "Protein.Q.Value",
            ],
        )
        if Protein_Qvalue:
            df = df[(df["Q.Value"] < 0.05) & (df["Protein.Q.Value"] < 0.05)]
        else:
            df = df[df["Q.Value"] < 0.05]
        # df = reduce_mem_usage(df, verbose=False)
        if df.shape[1] > 0:
            df.loc[:, "file_name"] = file.split("/")[-1]
    except Exception as e:
        print(
            f"\nWhile processing `{file}` file, \nThere occurred an error: \n {str(e)}\n"
        )
        return
    else:
        return df


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
    DIA_PATH = "/Users/ranpeng/Desktop/DIA-QC/data/likai-diann-raw/"
    # OUT_PATH = create_output_dir('results/QC_repeat90')
    # data = ThreadPool_load_file(inpath=DIA_PATH)
    data = joblib_load_file(inpath=DIA_PATH)
    print(data.shape)
    print(data.head())
    data = joblib_load_file_tqdm(inpath=DIA_PATH)
    print(data.shape)
    print(data.head())
    # # data = data.reset_index()
    # # print(data)
    # # data.to_feather(OUT_PATH + 'read_csv_joblib1.feather')
    # fullname_list = find_all_files_list(DIA_PATH, extension = '.tsv')
    # # print(fullname_list)
    # df = dd.read_csv(fullname_list)
    # print(df.shape)
    print(NUM_WORKERS)
