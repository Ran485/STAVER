#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File     : get_low_CI_peptides.py
@Time     : 2023/09/13 14:23:02
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
"""
# here put the import lib
import os
import pandas as pd
import numpy as np

from HighCIPeptides import HighCIPeptides
from staver.data import *
from utils import *


def preprocess_raw_peptides(data: pd.DataFrame, suffix="_F1_R1") -> pd.DataFrame:
    # Replace all the 0 values
    data = data.replace(0, np.nan)
    # data["file_name"] = data["file_name"].str.replace(
    #     r"DIA_\d*_HFX", "DIA_HFX", regex=True
    # )
    pattern = r"^.*(?=" + suffix + ")"
    data["file_name"] = data["file_name"].str.replace(
        pattern, "Reference_Benchmark_Dataset", regex=True
    )
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
    data = data.reset_index()
    # peptide_ident_single_file = data[data["count"] == 1].reset_index()
    # data = data[data["count"] > ount_cutoff_of_same_files].reset_index()
    # data = data[data["cv"] <= cv_thresh_peptide_of_same_files]
    data = index_transform(data)
    data = data.reset_index()

    return data


def get_low_CI_peptides():
    # load all reference peptide data
    data = pd.read_feather(
        "/Volumes/Disk3/STAVER-debug/preprocess-1/merge_all_peptide.feather"
    )
    data = preprocess_raw_peptides(data)
    # load high confidence peptide data
    high_CI_same_files = pd.read_csv(
        "/Volumes/Disk3/STAVER-debug/preprocess-1/peptides/high_confidence_peptide_of_same_files.csv"
    )
    high_CI_diff_files = pd.read_csv(
        "/Volumes/Disk3/STAVER-debug/preprocess-1/peptides/high_confidence_peptide_of_diff_files.csv"
    )

    high_CI = pd.concat([high_CI_same_files, high_CI_diff_files])

    data.set_index(["index", "file_name"], inplace=True)
    high_CI.set_index(["index", "file_name"], inplace=True)

    low_CI = data[~data.index.isin(high_CI.index)]
    low_CI.reset_index().to_csv(
        "/Volumes/Disk3/STAVER-debug/preprocess-1/low_CI_peptide-1.csv"
    )


if __name__ == "__main__":
    get_low_CI_peptides()
