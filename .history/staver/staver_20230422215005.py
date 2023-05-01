#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File     : staver_pipeline.py
@Time     : 2023/03/28 21:24:35
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2019-2020, DingLab-CHINA-SHNAGHAI
@Function : None
"""
# here put the import lib
import os
from HighCIPeptides import HighCIPeptides
from PeptideToProteinInference import PeptideToProteinInference
from merge_data import merge_data
from utils import create_output_dir, timer

# =======  只需要安装相关的包和配置以下信息  ======= #

# Setting the number of thresholds for computer operations
THRESHOD_NUMS = 16
# Whether to use standarde reference dataset
WITH_STANDARDE_DATA = True
# Setting the count cutoff of same files
COUNT_CUTOFF_SAME_FILES = 1
# Setting the count cutoff of different files
COUNT_CUTOFF_DIFF_FILES = 2
# Setting coefficient of variation threshold of the peptides
PEPTIDES_CV_THRESH = 0.3
# Setting the top high confidence interval peptides
TOP_PEPTIDES_NUMS = 15

# # DIA data directory
DIA_PATH = None  # "/Users/ranpeng/Desktop/DIA-QC/data/likai-diann-raw-20/"
# # DIA standarde reference directory
STANDARDE_DATASET_PATH = "/Users/ranpeng/Desktop/DIA-QC/data/likai-diann-raw/"
# The QC peptide data output directory
# create_output_dir('results/DIA_repeat201/peptides')
DIA_PEPTIDE_DATA_OUTPATH = None
# Standarde reference dataset
# REFERENCE_DATASET = "/Users/ranpeng/Desktop/DIA-QC/results/DIA_repeat20/peptides/\
#                      high_confidence_peptide_of_diff_files_top3_coefficient_top3.csv"
# The QC protein data output directory
# create_output_dir('results/DIA_repeat201/proteins')
DIA_PROTEIN_DATA_OUTPATH = None
# Setting coefficient of variation threshold for the proteins
PROTEINS_CV_THRESH = 0.3
# Setting the minimum threshold for a NUll peptide
NA_THRESHOLD = 0.3
# Set the suffix for folder specific identification
FILE_SUFFIX = "_F1_R1"


# Create an instance of the HighCIPeptides class
highCIPeptides = HighCIPeptides()
# Create an instance of the PeptideToProteinInference class
peptideToProteinInference = PeptideToProteinInference()


def preprocess_data(inpath: str, file_suffix) -> None:
    """Preprocess the data."""
    try:
        outpath = os.path.join(
            os.path.dirname(os.path.dirname(inpath)), "processed_data/"
        )
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        print(f"\nPreprocessed data stored in `{outpath}`\n")
        file_extend = f"{file_suffix}.csv"
        merge_data(
            inpath, usecols=["Gene_Symbol", "median"], filetype_extension=file_extend
        ).to_csv(f"{outpath}Protein_matrix.csv", index=False)

    except Exception as e:
        print(f"\nWhile merge`{inpath}`\n, There occurred an error:\n + {str(e)}\n")


@timer
def staver(
    number_threshods: int = THRESHOD_NUMS,
    with_standard_data: bool = True,
    dia_path: str = DIA_PATH,
    dia_pep_data_outpath: str = DIA_PEPTIDE_DATA_OUTPATH,
    standarde_dataset_path: str = STANDARDE_DATASET_PATH,
    file_suffix: str = FILE_SUFFIX,
    count_cutoff_same_files: int = COUNT_CUTOFF_SAME_FILES,
    count_cutoff_diff_files: int = COUNT_CUTOFF_DIFF_FILES,
    peptides_cv_thresh: float = PEPTIDES_CV_THRESH,
    # reference_dataset: str = REFERENCE_DATASET,
    dia_protein_data_outpath: str = DIA_PROTEIN_DATA_OUTPATH,
    proteins_cv_thresh: float = PROTEINS_CV_THRESH,
    top_peptides_nums: int = TOP_PEPTIDES_NUMS,
    na_threshold: float = NA_THRESHOLD,
) -> None:
    """Main function of STAVER algorithm operation.

    This function is the main function of STAVER algorithm operation.
    It is used to call the functions of each module to complete the
    operation of the STAVER algorithm.
    """
    # Save the log information to the log file
    # save_log(outpath=dia_protein_data_outpath)
    if with_standard_data:
        try:
            # Step 1: Get high confidence peptides
            highCIPeptides.num_workers = number_threshods
            highCIPeptides.reference_dataset_path = standarde_dataset_path
            highCIPeptides.count_cutoff_same_libs = count_cutoff_same_files
            highCIPeptides.count_cutoff_diff_libs = count_cutoff_diff_files
            highCIPeptides.cv_thresh_of_same_libs = peptides_cv_thresh
            highCIPeptides.cv_thresh_of_diff_libs = peptides_cv_thresh
            highCIPeptides.top_precursor_ions = top_peptides_nums
            highCIPeptides.outpath = dia_pep_data_outpath
            reference_dataset = highCIPeptides.main()

            # Step 2: Get high confidence proteins
            peptideToProteinInference.file_suffix = file_suffix
            peptideToProteinInference.input_data_path = dia_path
            peptideToProteinInference.reference_dataset = reference_dataset
            peptideToProteinInference.na_threshold = na_threshold
            peptideToProteinInference.cv_threshold = proteins_cv_thresh
            peptideToProteinInference.outpath = dia_protein_data_outpath
            peptideToProteinInference.parallel_main()

            # Step 3: Merge data
            preprocess_data(dia_protein_data_outpath, file_suffix)

        except Exception as e:
            error_message = (
                "\nWhile processing `"
                + DIA_PATH
                + " WITH_STANDARDE_DATA`\n, There occurred an error:\n"
                + str(e)
            )
            print(error_message)
    else:
        try:
            highCIPeptides.num_workers = number_threshods
            print("Load the default standard dataset for calibration!")
            # Step 2: Get high confidence proteins
            peptideToProteinInference.file_suffix = file_suffix
            peptideToProteinInference.input_data_path = dia_path
            peptideToProteinInference.outpath = dia_protein_data_outpath
            peptideToProteinInference.na_threshold = na_threshold
            peptideToProteinInference.cv_threshold = proteins_cv_thresh
            peptideToProteinInference.parallel_main()

            # Step 3: Merge data
            preprocess_data(dia_protein_data_outpath, file_suffix)

        except Exception as e:
            error_message = (
                "\nWhile processing `"
                + DIA_PATH
                + " WITHOUT_STANDARDE_DATA`\n, There occurred an error:\n"
                + str(e)
            )
            print(error_message)


if __name__ == "__main__":
    # The QC peptide data output directory
    DIA_PEPTIDE_DATA_OUTPATH = create_output_dir("results/DIA_rep20/peptides")
    # The QC protein data output directory
    DIA_PROTEIN_DATA_OUTPATH = create_output_dir("results/DIA_rep20/proteins")
    staver(
        with_standard_data=True,
        dia_pep_data_outpath=DIA_PEPTIDE_DATA_OUTPATH,
        dia_protein_data_outpath=DIA_PROTEIN_DATA_OUTPATH,
    )
