#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : HighCIPeptides.py
@Time     : 2023/05/02 00:07:09
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
'''
# here put the import lib


import pandas as pd
import numpy as np
import os
import gc
from utils import timer, status_info, reduce_mem_usage
from typing import Optional
from joblib import Parallel, delayed
from rich.progress import track
from read_file import *


"""
High Confidence Peptide Identification

This script identifies high confidence peptides from mass spectrometry data
by calculating the coefficient of variation (CV) of the peptide's peak area
and ranking the peptides based on their CV values. It then selects the top
precursor ions for further analysis.
"""


class HighCIPeptides:
    REFERENCE_DATASET_PATH: str = "./data/likai-diann-raw"
    COUNT_CUTOFF_OF_SAME_LIBS: int = 1
    COUNT_CUTOFF_OF_DIFF_LIBS: int = 2
    CV_THRESH_OF_SAME_LIBS: float = 0.3
    CV_THRESH_OF_DIFF_LIBS: float = 0.3
    TOP_PRECURSOR_IONS: int = 15
    OUTPATH = None  # create_output_dir('results/peptides')

    def __init__(
        self,
        num_workers: Optional[int] = NUM_WORKERS,
        reference_dataset_path: Optional[str] = REFERENCE_DATASET_PATH,
        count_cutoff_same_libs: Optional[int] = COUNT_CUTOFF_OF_SAME_LIBS,
        count_cutoff_diff_libs: Optional[int] = COUNT_CUTOFF_OF_DIFF_LIBS,
        cv_thresh_of_same_libs: Optional[float] = CV_THRESH_OF_SAME_LIBS,
        cv_thresh_of_diff_libs: Optional[float] = CV_THRESH_OF_DIFF_LIBS,
        top_precursor_ions: Optional[int] = TOP_PRECURSOR_IONS,
        outpath: Optional[str] = OUTPATH,
    ) -> None:
        self.num_workers = num_workers
        self.reference_dataset_path = reference_dataset_path
        self.count_cutoff_of_same_libs = count_cutoff_same_libs
        self.count_cutoff_of_diff_libs = count_cutoff_diff_libs
        self.cv_thresh_of_same_libs = cv_thresh_of_same_libs
        self.cv_thresh_of_diff_libs = cv_thresh_of_diff_libs
        self.top_precursor_ions = top_precursor_ions
        self.outpath = outpath
        """
        A class used to represent high confidence peptides.

        This script identifies high confidence peptides from mass spectrometry data
        by calculating the coefficient of variation (CV) of the peptide's peak area
        and ranking the peptides based on their CV values. It then selects the top
        precursor ions for further analysis.
        ...

        Parameters
            -------------
            reference_dataset_path : str
                The default path of the reference dataset.
            count_cutoff_of_same_libs : int
                The default count cutoff of the same libraries.
            count_cutoff_of_diff_libs : int
                The default count cutoff of different libraries.
            cv_thresh_of_same_libs : float
                The default CV threshold of the same libraries.
            cv_thresh_of_diff_libs : float
                The default CV threshold of different libraries.
            top_precursor_ions : int
                The default number of top precursor ions.
            outpath : str
                The default output path.

        Methods
        -------------
            load_data():
                Load the data from the input path.

            index_transform(df, convert_reverse=False):
                Transform the index of the given DataFrame.

            preprocess_filtered_data(data, ount_cutoff_of_same_files=None, 
                cv_thresh_peptide_of_same_files=None):
                Preprocess the filtered data.

            calculate_CV_for_same_files(data, count_cutoff_of_same_files=None, 
                cv_thresh_peptide_of_same_files=None):
                Calculate the coefficient of variation (CV) for the same files.

            calculate_CV_for_diff_files(df, ues_col='mean'):
                Calculate the coefficient of variation (CV) for different files.

            MaxSubsetMinCV(array):
                Find the subgroup of array with the lowest cv.

            high_qual_peptide_CV(df):
                Calculate the high quality peptide CV.

            find_subset(df):
                Find the subset in the DataFrame.

            applyParallel(dfGrouped, func):
                Apply the given function to each group in parallel.

            preprocess_raw_data(data, replace_0=True):
                Preprocess the raw data.

            identify_high_confidence_peptides(df):
                Generate a dictionary of high confidence peptides.

            get_peptide_coefficient(df):
                Get the peptide coefficient.

            get_top_high_ci_peptides(df, top_precursor_ion_nums=None):
                Get the top high confidence peptides.
        """

    @timer
    def load_data(self) -> pd.DataFrame:
        """
        Load the data from the specified input data path.

        Returns
        -------
        pd.DataFrame
            The loaded data as a pandas DataFrame.
        """
        if os.path.exists(self.outpath + "merge_all_peptide.feather"):
            data = pd.read_feather(self.outpath + "merge_all_peptide.feather")
        else:
            print(self.outpath)
            data = joblib_load_file(
                inpath=self.reference_dataset_path,
                extension=".tsv",
                num_workers=self.num_workers,
            )
            data = reduce_mem_usage(data, verbose=False)
            print(
                f"\nThe merged raw data save to the follow directory: \n`{self.outpath}merge_all_peptide.feather`"
            )
            data.reset_index().to_feather(
                self.outpath + "merge_all_peptide.feather", compression="lz4"
            )
        return data

    @staticmethod
    def index_transform(df, convert_reverse=False) -> pd.DataFrame:
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
            )
            df.set_index("index", inplace=True)
        return df

    def preprocess_filtered_data(
        self,
        data: pd.DataFrame,
        ount_cutoff_of_same_files=None,
        cv_thresh_peptide_of_same_files=None,
    ) -> pd.DataFrame:
        # Replace all the 0 values
        data = data.replace(0, np.nan)
        peptides = data.groupby(
            ["Genes", "Stripped.Sequence", "Precursor.Id", "Modified.Sequence"]
        )["Precursor.Normalised"].agg(["mean", "std", "count"])
        peptides["cv"] = peptides["std"] / peptides["mean"]
        # peptides.to_csv(self.outpath + "peptide_cv_20221207.csv")
        peptides_num = peptides.shape[0]
        data = data.groupby(
            [
                "file_name",
                "Genes",
                "Stripped.Sequence",
                "Precursor.Id",
                "Modified.Sequence",
            ]
        )["Precursor.Normalised"].agg(["mean", "std", "count"])
        data["cv"] = data["std"] / data["mean"]
        peptide_ident_single_file = data[data["count"] == 1].reset_index()
        data = data[data["count"] > ount_cutoff_of_same_files].reset_index()
        data = data[data["cv"] <= cv_thresh_peptide_of_same_files]
        data = self.index_transform(data)
        peptide_ident_single_file = self.index_transform(peptide_ident_single_file)
        return data, peptide_ident_single_file, peptides_num

    def calculate_CV_for_same_files(
        self,
        data: pd.DataFrame,
        count_cutoff_of_same_files=None,
        cv_thresh_peptide_of_same_files=None,
    ) -> pd.DataFrame:
        # Calculate the CV for the same files
        data, peptide_ident_single_file, peptides_num = self.preprocess_filtered_data(
            data, count_cutoff_of_same_files, cv_thresh_peptide_of_same_files
        )
        filtered_highCI_peptides_num = (
            data.groupby(data.index)["mean"].agg(["mean", "count"]).shape[0]
        )
        filtered_onefile_peptides_num = (
            peptide_ident_single_file.groupby(peptide_ident_single_file.index)["mean"]
            .agg(["mean", "count"])
            .shape[0]
        )
        # The number of peptides that are identified in only one file
        filtered_peptides_num = (
            filtered_highCI_peptides_num + filtered_onefile_peptides_num
        )
        high_confidence_peptides_ratio = filtered_highCI_peptides_num / peptides_num
        onefile_peptides_ratio = filtered_onefile_peptides_num / peptides_num
        print(f"The high confidence peptides number is: {filtered_highCI_peptides_num}")
        print(
            f"The onefile identified peptides number is: {filtered_onefile_peptides_num}"
        )
        print(
            f"The high confidence peptides ratio is: {high_confidence_peptides_ratio}%"
        )
        print(f"The onefile peptides ratio is: {onefile_peptides_ratio}%")
        peptide_loss_rate = 1 - (filtered_highCI_peptides_num / peptides_num)
        print(
            f"The high CI peptides of the same file filter loss rate is: {peptide_loss_rate * 100}%\n"
        )
        data_processed = data[["file_name", "mean", "cv", "count"]]
        del data
        gc.collect()
        return data_processed, peptide_ident_single_file, filtered_peptides_num

    @staticmethod
    def calculate_CV_for_diff_files(df, ues_col="mean") -> float:
        mean = df[ues_col].mean()
        std = df[ues_col].std()
        cv = std / mean
        return cv, mean

    def MaxSubsetMinCV(self, array: pd.DataFrame) -> pd.DataFrame:
        """Find the subgroup of array with the lowest cv."""
        if len(array) == 0:
            return None
        elif len(array) == 1:
            return array
        else:
            # array.sort(reverse=True)
            array = array.sort_values("mean", ascending=False)
            curr_cv, _ = self.calculate_CV_for_diff_files(array)
            minimum_cv = curr_cv
            # Loop through the data frame
            for row in range(len(array)):
                # If the number of rows is less than or equal to count_cutoff, break the loop
                if (
                    len(array) <= self.count_cutoff_of_diff_libs
                    or minimum_cv <= self.cv_thresh_of_diff_libs
                ):
                    break
                # Pop the last row from the data frame
                array = array.iloc[:-1]
                curr_cv, _ = self.calculate_CV_for_diff_files(array)
                if curr_cv < minimum_cv:
                    minimum_cv = curr_cv
            # print(f'\nThe minimum cv is: {minimum_cv}')
            # print(f'\nThe minimum cv list is: {array}')
            if minimum_cv <= self.cv_thresh_of_diff_libs:
                return array

    def high_qual_peptide_CV(self, df) -> pd.DataFrame:
        _df = df.groupby(df.index)["mean"].agg(["mean", "std", "count"])
        _df["cv"] = _df["std"] / _df["mean"]
        return _df

    def find_subset(self, df) -> pd.DataFrame:
        if len(df) > 1:
            return df

    @staticmethod
    def applyParallel(dfGrouped, func):
        """
        Apply a given function to the data in parallel.

        Parameters
        ----------
        dfGrouped : pd.DataFrame
            The grouped input data.
        func : function
            The function to be applied in parallel.

        Returns
        -------
        pd.DataFrame
            The DataFrame with the applied function.
        """
        res = Parallel(n_jobs=NUM_WORKERS)(
            delayed(func)(group)
            for name, group in track(dfGrouped, description="Parallel processing...")
        )
        return pd.concat(res)

    @staticmethod
    def preprocess_raw_data(data, replace_0=True):
        data["file_name"] = data["file_name"].str.replace(
            r"DIA_\d*_HFX", "DIA_HFX", regex=True
        )
        # data['file_name'] = data['file_name'].apply(str_replace)
        # tmp.loc[:, 'file_name'] = tmp.loc[:, 'file_name'].str.split(pat=r'F1_R\d_', expand=True)[1]
        # tmp['file_name'] = tmp['file_name'].replace(r'min_\w*\d*_2per_\d*_F', 'min_2per_F', regex=True)
        if replace_0:
            data = data.replace(0.0, np.nan)
        data = data[data["Genes"].notnull()]
        processed_data = data[~data["Genes"].str.startswith("rev")]
        return processed_data

    @timer
    def identify_high_confidence_peptides(self, df) -> pd.DataFrame:
        """Generate a dictionary of high confidence peptides.

        Parameters
        ----------
        df : pd.DataFrame
            The input data as a pandas DataFrame.

        Returns
        -------
            high_confidence_peptide_dict: dict
                A dictionary of high confidence peptides.
        """
        # Process the input file
        df = self.preprocess_raw_data(df)
        # Drop the NAN values
        df = df.dropna(subset=["Precursor.Normalised"])
        # Calculate the coefficient of variation of the data
        df, peptide_ident_single_file, peptides_num = self.calculate_CV_for_same_files(
            df,
            count_cutoff_of_same_files=self.count_cutoff_of_same_libs,
            cv_thresh_peptide_of_same_files=self.cv_thresh_of_same_libs,
        )
        # Export the high_confidence_peptide of same files
        df.to_csv(self.outpath + "high_confidence_peptide_of_same_files.csv")
        peptide_ident_single_file.to_csv(
            self.outpath + "peptide_of_onefile_identify.csv"
        )
        high_confidence_peptides = self.applyParallel(
            df.groupby(df.index)[["file_name", "mean", "cv", "count"]],
            self.MaxSubsetMinCV,
        )
        # Discard those that have been identified only once.
        peptide_ident_single_file = self.applyParallel(
            peptide_ident_single_file.groupby(peptide_ident_single_file.index)[
                ["file_name", "mean", "cv", "count"]
            ],
            self.MaxSubsetMinCV,
        )
        high_confidence_peptides_onefile = self.applyParallel(
            peptide_ident_single_file.groupby(peptide_ident_single_file.index)[
                ["file_name", "mean", "cv", "count"]
            ],
            self.find_subset,
        )
        # Export the high_confidence_peptide of different files
        high_confidence_peptides = pd.concat(
            [high_confidence_peptides, peptide_ident_single_file]
        )
        high_confidence_peptides.to_csv(
            self.outpath + "high_confidence_peptide_of_diff_files.csv"
        )
        del df
        gc.collect()
        high_confidence_peptides_cv = self.high_qual_peptide_CV(
            high_confidence_peptides
        )
        # Export the high_confidence_peptide cv of different files
        high_confidence_peptides_cv.to_csv(
            self.outpath + "high_confidence_peptide_of_diff_files_cv.csv"
        )
        peptide_onefile_idendify_of_diff_files = peptide_ident_single_file[
            ~peptide_ident_single_file.index.isin(
                high_confidence_peptides_onefile.index
            )
        ]
        # Export the peptide_onefile_idendify_of_diff_files of different files
        peptide_onefile_idendify_of_diff_files.to_csv(
            self.outpath + "peptide_of_onefile_identify_of_diff_files.csv"
        )
        high_CI_peptides_num = high_confidence_peptides_cv.shape[0]
        peptide_loss_rate = 1 - (high_CI_peptides_num / peptides_num)
        print(f"\nThe high CI peptides number is: {high_CI_peptides_num}")
        print(
            f"The number of one idedified peptides with no confidence is: {peptide_onefile_idendify_of_diff_files.shape[0]}"
        )
        print(
            f"The high CI peptides of the different files filter loss rate is: {peptide_loss_rate * 100}%"
        )
        high_CI_ratio = high_CI_peptides_num / (
            peptide_onefile_idendify_of_diff_files.shape[0] + high_CI_peptides_num
        )
        print(f"The ratio of high confidence peptides is: {high_CI_ratio}%\n")

        return high_confidence_peptides, high_confidence_peptides_cv

    def get_peptide_coefficient(self, df) -> pd.DataFrame:
        df["Coefficient"] = 1 / (df["mean"] / (df["mean"].sum() / len(df)))
        df["mean1"] = df["mean"] * df["Coefficient"]
        return df

    @timer
    def get_top_high_ci_peptides(self, df, top_precursor_ion_nums=None) -> pd.DataFrame:
        top_files = self.applyParallel(
            df.groupby(df.index)[["file_name", "mean", "cv", "count"]],
            lambda x: x.sort_values(["count", "cv"], ascending=[False, True]).head(
                top_precursor_ion_nums
            ),
        )
        top_files.to_csv(
            self.outpath + "high_confidence_peptide_of_diff_files_top3.csv"
        )
        top_files = top_files.reset_index()
        top_files["protein_index"] = (
            top_files["index"]
            .str.split("_", expand=True)[0]
            .str.split(";", expand=True)[0]
        )
        top_files = top_files.set_index("index")
        top_peptides = self.applyParallel(
            top_files.groupby("protein_index")[["file_name", "mean", "cv", "count"]],
            lambda x: x.sort_values(["count", "cv"], ascending=[False, True]).head(
                top_precursor_ion_nums
            ),
        )
        top_peptides = self.applyParallel(
            top_peptides.groupby(top_peptides.index)[
                ["file_name", "mean", "cv", "count"]
            ],
            self.get_peptide_coefficient,
        )
        top_peptides.to_csv(
            self.outpath + f"high_confidence_peptides_top{top_precursor_ion_nums}.csv"
        )
        reference_path = os.path.join(
            self.outpath, f"high_confidence_peptides_top{top_precursor_ion_nums}.csv"
        )
        return reference_path

    def main(self):
        """
        The main function of the HighCIPeptides class.

        Loads the mass spectrometry data, identifies high confidence peptides,
        and selects the top high confidence peptides based on their CV values.
        """
        data = self.load_data()
        high_confidence_peptides, _ = self.identify_high_confidence_peptides(data)
        reference_path = self.get_top_high_ci_peptides(
            high_confidence_peptides, self.top_precursor_ions
        )
        del data
        gc.collect()
        return reference_path


if __name__ == "__main__":
    # Initialize the defined class
    staver = HighCIPeptides()
    staver.main()
