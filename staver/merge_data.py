#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File     : merge_data.py
@Time     : 2022/11/10 19:09:39
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
"""
# here put the import lib

import pandas as pd
import numpy as np
import os
from rich.console import Console


def find_all_files(input_folder, extension=".csv"):
    """
    Finds all files in a given folder with a specific extension.

    Args:
        input_folder (str): Directory to search for files.
        extension (str): File extension to look for. Defaults to ".csv".

    Returns:
        list: A list of full file paths matching the given extension.
    """
    fullname = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(extension):
                path = os.path.join(root, file)
                fullname.append(path)
    return fullname


def read_file(file):
    """
    Reads a file based on its extension and returns its content as a DataFrame.

    Args:
        file (str): Path to the file to be read.

    Returns:
        pd.DataFrame: DataFrame containing the file's content.
    """
    if file.endswith(".csv"):
        return pd.read_csv(file)
    elif file.endswith(".txt"):
        return pd.read_table(file)
    elif file.endswith(".xlsx"):
        return pd.read_excel(file)
    else:
        raise ValueError(f"Unsupported file extension in file {file}")


def merge_data(path, filetype_extension=".txt", usecols=None):
    """
    Merges data from multiple files in a given directory into a single DataFrame.

    Args:
        path (str): Path to the directory containing the files.
        filetype_extension (str): Extension of files to merge. Defaults to '.txt'.
        usecols (list): Columns to use from the files. Defaults to ['Symbol', 'Area'].

    Returns:
        pd.DataFrame: DataFrame containing merged data from all files.
    """
    if usecols is None:
        usecols = ["Symbol", "Area"]
    files = find_all_files(path, extension=filetype_extension)
    console = Console()
    merged_df = pd.DataFrame()

    for file in files:
        try:
            file_data = read_file(file)
            file_data = file_data[usecols]
            file_data.rename(
                columns={file_data.columns[1]: file.split("/")[-1].split(".")[0]},
                inplace=True,
            )
            merged_df = pd.merge(
                merged_df,
                file_data,
                on=usecols[0],
                how="outer",
                suffixes=(False, False),
            )
            console.print(f"[green]Successfully merged file:[/green] {file}")
        except Exception as e:
            console.print(f"[red]Failed to merge file:[/red] {file}\nError: {e}")

    return merged_df


def deduplicate(data):
    """
    Deduplicates rows in the DataFrame and aggregates data.

    Args:
        data (pd.DataFrame): The DataFrame to deduplicate.

    Returns:
        pd.DataFrame: Deduplicated DataFrame.
    """
    console = Console()
    console.print("Before drop duplicate:", data.shape)
    console.print("The total duplicate number is:", data.duplicated().sum())
    data.drop_duplicates(inplace=True)
    data = data.groupby(data.index).agg("sum").replace(0, np.nan)
    console.print("After drop duplicate:", data.shape)
    return data


if __name__ == "__main__":
    path = r"/path/to/data/"
    outpath = r"/path/to/output/"
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    df = merge_data(path, filetype_extension=".csv", usecols=["Symbol", "Area"])
    df = deduplicate(df)
    df.to_csv(f"{outpath}merged_matrix_rawdata.csv", index=False)
