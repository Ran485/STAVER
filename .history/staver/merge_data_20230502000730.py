

import pandas as pd
import numpy as np
import os
from rich.console import Console

console = Console()


def find_all_files(input_folder, extension=".csv"):
    fullname = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(extension):
                path = os.path.join(root, file)
                fullname.append(path)
    return fullname


def read_files(file):
    if file.endswith(".csv"):
        file = pd.read_csv(file)
    elif file.endswith(".txt"):
        file = pd.read_table(file)
    elif file.endswith(".xlsx"):
        file = pd.read_excel(file)
    return file


def merge_data(
    path=None, filetype_extension=".txt", usecols=None, Unique_Peptide_Num=1
):
    """merge the download txt files from firmiana.

    Args:
    -------------
        path (str, optional):
                        [description]. Defaults to None.
        filetype_extension (str, optional):
                        [description]. Defaults to '.txt'.
        usecols (list, optional):
                        [description]. Defaults to ['Symbol', 'Area'].
        Unique_Peptide_Num (int, optional):
                        [description]. Defaults to 1.

    Returns:
    -------------
        dataframe: the merged csv dataframe of all the files in the path.
    """
    if usecols is None:
        usecols = ["Symbol", "Area"]
    files = find_all_files(path, extension=filetype_extension)
    filenames = [file.split(".")[0].split("/")[-1] for file in files]
    # print(filenames)
    merge_df = pd.DataFrame([], columns=[usecols[0]])
    for file, filename in zip(
        files, filenames
    ):  # , description = f"[bold][orange]Merging Data..."):
        try:
            # print(f"Start merging {file} samples!")
            file = read_files(file)
            # Filter the unique peptide number
            # file = file[file['Unique Peptide Num'] >= Unique_Peptide_Num]
            file = file[usecols]
            file.rename(columns={file.columns[1]: filename}, inplace=True)
            merge_df = pd.merge(merge_df, file, on=usecols[0], how="outer")
            # merge_data = merge_data.sort_values('Filename')
            console.print(f"[green]Successfully merged sample:[/green] {filename}!")
        except Exception as e:
            console.print(f"\n[bold red]Failed merge file:[/bold red] {filename}")
            print(f"Error in file: {e}\n")

    return merge_df


def deduplicate(data):
    # Drop duplicate rows
    print("Before drop duplicate:", data.shape)
    print("The toatal duplicate number is:", data.duplicated().sum())
    data.drop_duplicates(inplace=True)
    # Aggregate data by summing values for the same index
    data = data.groupby(data.index).agg("sum").replace(0, np.nan)
    print("After drop duplicate:", data.shape)
    return data


if __name__ == "__main__":
    path = r"/Volumes/T7_Shield/staver/results/DIA_RM/pep_cv03_rawdata_20230423_ge2/"
    # outpath = create_output_dir('results/merge_data')
    outpath = r"/Volumes/T7_Shield/staver/results/DIA_RM/"
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    df = merge_data(path, usecols=["Symbol", "Area"], filetype_extension=".csv")
    # df = deduplicate(df)
    df.to_csv(f"{outpath}Protein_matrix_rawdata.csv", index=False)
