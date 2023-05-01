#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : staver_pipeline.py
@Time     : 2023/05/02 00:08:46
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
'''
# here put the import lib


# Combine modules HighCIPeptides and PeptideToProteinInference into
# a list of parameters that can be called from the command line

import argparse
import os
from staver import staver
from utils import save_log


class ValidateParser:
    """Validate the input arguments.

    This class is used to validate the input arguments.

    Attributes:
        None

    Methods:
        valid_input_path: Validate the input path.
        valid_output_path: Validate the output path.
        valid_int: Validate the input int value.
        valid_float: Validate the input float value.
    """

    @staticmethod
    def valid_input_path(path):
        # Check if the path exists and is a directory
        if os.path.exists(path) and os.path.isdir(path):
            return path
        else:
            raise argparse.ArgumentTypeError(
                f"The path: '{path}' does not exist.\n \
                Please check the path and provide an valid input path."
            )

    @staticmethod
    def valid_output_path(path):
        # Check if the path exists and is a directory
        if not os.path.exists(path):
            os.makedirs(path)
            print(f"Sucessfully created the output directory: {path}")
            return path
        elif os.path.exists(path) and os.path.isdir(path):
            return path
        else:
            raise argparse.ArgumentTypeError(
                f"The path: '{path}' is not a directory.\n \
                Please check the path and provide an valid output directory."
            )

    @staticmethod
    def valid_int(value):
        try:
            int_value = int(value)
        except ValueError as ve:
            raise argparse.ArgumentTypeError(
                f"Invalid int value: '{value}', must be an integer."
            ) from ve

        if int_value < 0:
            raise argparse.ArgumentTypeError(
                f"Invalid int value: '{value}', \
                must be an integer greater than 0."
            )
        return int_value

    @staticmethod
    def valid_float(value):
        try:
            float_value = float(value)
        except ValueError as ve:
            raise argparse.ArgumentTypeError(
                f"Invalid float value: '{value}', must be a float between 0-1."
            ) from ve

        if float_value < 0 or float_value > 1:
            raise argparse.ArgumentTypeError(
                f"Invalid float value: '{value}', must be a float between 0-1."
            )
        return float_value


def print_arguments(args):
    print("All parsed arguments:")
    # Print all the arguments in the dictionary
    for arg in vars(args):
        print(f"{arg}: {getattr(args, arg)}")


def export_arguments(args, output_path):
    with open(os.path.join(output_path, "staver_arguments.txt"), "w") as f:
        f.write("All parsed arguments:  \r ")
        # Print all the arguments in the dictionary
        for arg in vars(args):
            f.write(f"{arg}: {getattr(args, arg)}  \r ")


def staver_pipeline():
    """Main function of STAVER algorithm operation.

    This function is the main function of STAVER algorithm operation.
    It is used to call the functions of each module to complete the
    operation of the STAVER algorithm.
    """

    parser = argparse.ArgumentParser(
        description="STAVER: A Standardized Dataset-Based Algorithm for \
                    Efficient Variation Reduction in Large-Scale DIA MS Data"
    )
    parser.add_argument(
        "-n",
        "--thread_numbers",
        required=True,
        dest="number_threshods",
        type=ValidateParser.valid_int,
        help="The number of thresholds for computer operations",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        dest="dia_path",
        type=ValidateParser.valid_input_path,
        help="The DIA input data path",
    )
    parser.add_argument(
        "-ref",
        "--reference_dataset_path",
        required=False,
        type=ValidateParser.valid_input_path,
        dest="reference_standard_dataset",
        help="The DIA standarde reference directory",
    )
    parser.add_argument(
        "-o",
        "--output_peptide",
        required=True,
        type=ValidateParser.valid_output_path,
        dest="dia_pep_data_outpath",
        help="The processed DIA proteomics of peptide data output path",
    )
    parser.add_argument(
        "-op",
        "--output_protein",
        required=True,
        type=ValidateParser.valid_output_path,
        dest="dia_protein_data_outpath",
        help="The processed DIA proteomics protein data output path",
    )
    parser.add_argument(
        "-c",
        "--count_cutoff_same_libs",
        required=False,
        type=ValidateParser.valid_int,
        default=1,
        help="Setting the count cutoff of same files (default: 1)",
    )
    parser.add_argument(
        "-d",
        "--count_cutoff_diff_libs",
        required=False,
        type=ValidateParser.valid_int,
        default=2,
        help="Setting the count cutoff of different files (default: 2)",
    )
    parser.add_argument(
        "-pep_cv",
        "--peptides_cv_thresh",
        required=False,
        type=ValidateParser.valid_float,
        default=0.3,
        help="Setting coefficient of variation threshold \
            for the peptides (default: 0.3)",
    )
    parser.add_argument(
        "-pro_cv",
        "--proteins_cv_thresh",
        required=False,
        type=ValidateParser.valid_float,
        default=0.3,
        help="Setting coefficient of variation threshold \
            for the proteins (default: 0.3)",
    )
    parser.add_argument(
        "-na_thresh",
        "--na_threshold",
        required=False,
        type=ValidateParser.valid_float,
        default=0.3,
        help="Setting the minimum threshold for \
                NUll peptides (default: 0.3)",
    )
    parser.add_argument(
        "-top",
        "--top_precursor_ions",
        required=False,
        type=ValidateParser.valid_int,
        default=10,
        help="Setting the top high confidence interval precursor ions",
    )
    parser.add_argument(
        "-suffix",
        "--file_suffix",
        required=False,
        type=str,
        default="_F1_R1",
        help="Set the suffix for folder specific identification",
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0.0",
    )

    # Print and export the parsed arguments to a dictionary
    args = parser.parse_args()
    # Save the log information to the log file
    # save_log(outpath=args.dia_protein_data_outpath)
    print_arguments(args)
    export_arguments(args, args.dia_pep_data_outpath)

    if args.reference_standard_dataset is not None:
        # Step 1: Get high confidence peptides
        staver(
            number_threshods=args.number_threshods,
            with_standard_data=True,
            dia_path=args.dia_path,
            dia_pep_data_outpath=args.dia_pep_data_outpath,
            standarde_dataset_path=args.reference_standard_dataset,
            file_suffix=args.file_suffix,
            count_cutoff_same_files=args.count_cutoff_same_libs,
            count_cutoff_diff_files=args.count_cutoff_diff_libs,
            peptides_cv_thresh=args.peptides_cv_thresh,
            proteins_cv_thresh=args.proteins_cv_thresh,
            dia_protein_data_outpath=args.dia_protein_data_outpath,
            top_peptides_nums=args.top_precursor_ions,
            na_threshold=args.na_threshold,
        )

    else:
        print("Load the default standard dataset for calibration!")
        print(args.dia_path)
        staver(
            number_threshods=args.number_threshods,
            with_standard_data=False,
            dia_path=args.dia_path,
            dia_pep_data_outpath=args.dia_pep_data_outpath,
            file_suffix=args.file_suffix,
            count_cutoff_same_files=args.count_cutoff_same_libs,
            count_cutoff_diff_files=args.count_cutoff_diff_libs,
            peptides_cv_thresh=args.peptides_cv_thresh,
            proteins_cv_thresh=args.proteins_cv_thresh,
            dia_protein_data_outpath=args.dia_protein_data_outpath,
            top_peptides_nums=args.top_precursor_ions,
            na_threshold=args.na_threshold,
        )


if __name__ == "__main__":
    staver_pipeline()
