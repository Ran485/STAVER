#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File     : get_high_confidence_proteins.py
@Time     : 2022/09/28 14:56:50
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2019-2020, DingLab-CHINA-SHNAGHAI
@Function : None
"""
# here put the import lib
import pandas as pd
import os
from typing import Optional
from HighCIPeptides import *
from utils import *
from data import *
from batch_rename import batch_renaming
from MaxLFQ import *


"""
The Peptide To Protein Inference

This Python script defines a class called PeptideToProteinInference, which 
is designed to infer proteins from peptides with high confidence. The class
takes in several optional arguments such as file_suffix, input_data_path, 
reference_dataset, na_threshold, cv_threshold, and outpath. The class contains
several methods to perform different tasks in the inference process.
"""


class PeptideToProteinInference:
    FILE_SUFFIX = "_F1_R2"
    INPUT_DATA_PATH: str = "./data/likai-diann-raw-20/"
    REFERENCE_DATASET: str = "/Volumes/T7_Shield/staver/results/DIA_repeat20/peptides/high_confidence_peptides_top3.csv"
    NA_THRESHOLD: float = 0.3
    CV_THRESHOLD: float = 0.1
    OUTPATH: str = None  # create_output_dir('results/proteins')

    def __init__(
        self,
        num_workers: Optional[int] = NUM_WORKERS,
        file_suffix: Optional[str] = FILE_SUFFIX,
        input_data_path: Optional[str] = INPUT_DATA_PATH,
        reference_dataset: Optional[str] = REFERENCE_DATASET,
        na_threshold: Optional[float] = NA_THRESHOLD,
        cv_threshold: Optional[float] = CV_THRESHOLD,
        outpath: Optional[str] = OUTPATH,
    ) -> pd.DataFrame:
        self.num_workers = num_workers
        self.file_suffix = file_suffix
        self.input_data_path = input_data_path
        self.reference_dataset = reference_dataset
        self.na_threshold = na_threshold
        self.cv_threshold = cv_threshold
        self.outpath = outpath
        """
        The Peptide To Protein Inference

        This Python script defines a class called PeptideToProteinInference, which 
        is designed to infer proteins from peptides with high confidence. The class
        takes in several optional arguments such as file_suffix, input_data_path, 
        reference_dataset, na_threshold, cv_threshold, and outpath. The class contains
        several methods to perform different tasks in the inference process.
        
        Parameters:
            -----------
            file_suffix : str, optional
                The suffix of the input data file. (default is "_F1_R1")
            input_data_path : str, optional
                The path of the input data. (default is "/Volumes/T7_Shield/staver/data/likai-diann-raw-20/")
            reference_dataset : str, optional
                The path of the reference dataset. (default is "/Volumes/T7_Shield/staver/results/peptides/high_confidence_peptide_of_diff_files_top3_coefficient_top3.csv")
            na_threshold : float, optional
                The threshold of the missing value. (default is 0.3)
            cv_threshold : float, optional
                The threshold of the coefficient of variation. (default is 0.1)
            outpath : str, optional
                The path of the output data. (default is "/Volumes/T7_Shield/staver/results/proteins/")

        Returns:
            -----------
            CSV files contains the inferred proteins.
            """

    def find_all_foder(self) -> list:
        #  find all directory in a directory
        folders = []
        for root, dirs, _ in os.walk(self.input_data_path):
            folders.extend(os.path.join(root, dir) for dir in dirs)
        return folders

    def fast_scandir(self, dirname) -> list:
        subfolders = [f.path for f in os.scandir(dirname) if f.is_dir()]
        for dirname in list(subfolders):
            subfolders.extend(self.fast_scandir(dirname))
        return subfolders

    def load_high_ci_peptides(self) -> pd.DataFrame:
        print("Loading high confidence peptides form ")
        print(self.reference_dataset)
        df = pd.read_csv(self.reference_dataset)
        df["file_name"] = df["file_name"].str.split("_F1_R1_", expand=True)[1]
        df["index"] = df["index"] + "__" + df["file_name"]
        df = df[["index", "Coefficient"]]
        return df

    def caculate_CV_for_filtered_peptides(self, df) -> pd.DataFrame:
        df = (
            df.groupby("index")["Precursor.Normalised"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )
        df["cv"] = df["std"] / df["mean"]
        return df

    def get_high_CI_peptides(self, df, cv_thresh=None) -> pd.DataFrame:
        df = df[(df["cv"] < cv_thresh) | (df["count"] == 1)]
        return df

    def FilterMinCV(self, array) -> pd.DataFrame:
        """Find the subgroup of array with the lowest cv."""
        if len(array) == 0:
            return None
        elif len(array) == 1:
            return array
        else:
            # array.sort(reverse=True)
            array = array.sort_values(
                by=["mean", "count", "cv"], ascending=[False, False, True]
            )
            array_min_cv = array[array["cv"] < self.cv_threshold]
            if len(array_min_cv) >= 1:
                return array_min_cv
            row_num = round(len(array) * 0.2, 0)
            array_min_cv = array[: int(row_num)]
            if array_min_cv["cv"].mean() < 0.8:
                return array[: int(row_num)]

    def get_high_confidence_proteins(self, data, high_ci_peptides) -> pd.DataFrame:
        # Step 1: Make a copy of the input data to work with
        d = data.copy()

        # Step 2: Method chaining to avoid intermediate variables
        suffix = self.file_suffix + "_"
        peptides = (
            d.pipe(reduce_mem_usage, verbose=False)
            .pipe(index_transform)
            .reset_index()
            .assign(
                file_name=lambda x: x["file_name"].str.split(suffix, expand=True)[1],
                index=lambda x: x["index"] + "__" + x["file_name"],
            )
            .merge(high_ci_peptides, on="index", how="inner")
        )

        data_ = (
            peptides.copy()
            .assign(index=lambda x: x["index"].str.rsplit("__", 1, expand=True)[0])
            .assign(
                Precursor_Normalized=lambda x: x["Precursor.Normalised"]
                * x["Coefficient"]
            )
            .pipe(self.caculate_CV_for_filtered_peptides)
            .assign(
                Gene_Symbol=lambda x: x["index"]
                .str.split("_", expand=True)[0]
                .str.split(";", expand=True)[0]
            )
            .set_index("Gene_Symbol")  # in-place indexing is slower than assignment
        )
        # !Step 3: Use groupby-apply instead of applyParallel
        # data_ = data_.groupby(data_.index).apply(self.FilterMinCV).reset_index()
        data_ = applyParallel(
            data_.groupby(data_.index), self.FilterMinCV
        ).reset_index()

        # Step 4: Remove unnecessary operations
        data_ = data_.groupby("Gene_Symbol")["mean"].agg(
            ["sum", "count", "mean", "median"]
        )
        data_.rename(
            {"sum": "Abundance", "count": "Unique_peptides"}, axis=1, inplace=True
        )
        return data_, peptides

    def main(self, subdirectory, high_ci_peptides, outpath=None) -> pd.DataFrame:
        if os.path.basename(subdirectory).endswith(self.file_suffix):
            # try:
            data = joblib_load_file(subdirectory)
            data, peptides = self.get_high_confidence_proteins(data, high_ci_peptides)
            data.to_csv(
                os.path.join(outpath, os.path.basename(subdirectory) + ".csv"),
                index=True,
            )
            peptides.to_csv(
                os.path.join(outpath, os.path.basename(subdirectory) + "_peptides.csv"),
                index=False,
            )
            # except Exception as e:
            #     print(
            #         f"\nwhile processing foder: `{os.path.basename(subdirectory)}` , there occurred an error: \n + {str(e)}"
            #     )
            #     return None

    # @status_info()
    def parallel_main(self) -> None:
        high_ci_peptides = self.load_high_ci_peptides()
        folders = self.fast_scandir(self.input_data_path)
        self.check_file_suffix(folders)

        folders = self.fast_scandir(self.input_data_path)
        # Increase worker timeout time to avoid memory leaks
        Parallel(n_jobs=self.num_workers - 4, timeout=3600)(
            delayed(self.main)(subdirectory, high_ci_peptides, self.outpath)
            for subdirectory in folders
            if os.path.basename(subdirectory).endswith(self.file_suffix)
        )

    def check_file_suffix(self, folders):
        if self.file_suffix not in folders[0] or "__" in folders[1]:
            batch_renaming(self.input_data_path, self.file_suffix)


if __name__ == "__main__":
    # Initialize the defined class
    FILE_SUFFIX = "_F1_R2"
    PepProInfer = PeptideToProteinInference()
    PepProInfer.file_suffix = "_F1_R2"
    PepProInfer.input_data_path = "/Volumes/Samsung_T5/DIA-QC/data/file_transfer/demo/"
    PepProInfer.outpath = create_output_dir("results/demo/")
    # PepProInfer.file_suffix = "_F1_R2"
    PepProInfer.parallel_main()
