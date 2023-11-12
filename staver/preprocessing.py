#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@File     : preprocessing.py
@Time     : 2022/11/10 10:26:18
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
    """
    Reduce the memory usage of a pandas DataFrame by downcasting numeric types.

    Iterates over each column and determines if the numeric data can be
    downcasted to a more memory-efficient data type, such as int or float
    with lower precision.

    Args:
        df (pd.DataFrame): The DataFrame whose memory usage is to be reduced.
        verbose (bool, optional): Whether to print the memory reduction information. Defaults to True.

    Returns:
        pd.DataFrame: DataFrame with optimized memory usage.

    Examples:
        >>> data = {'int_column': [1, 2, 3], 'float_column': [0.1, 0.2, 0.3]}
        >>> df = pd.DataFrame(data)
        >>> original_memory = df.memory_usage().sum() / 1024**2
        >>> print(f'Original memory usage: {original_memory:.2f} MB')
        >>> df_reduced = reduce_mem_usage(df)
        >>> reduced_memory = df_reduced.memory_usage().sum() / 1024**2
        >>> print(f'Reduced memory usage: {reduced_memory:.2f} MB')
    """
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
    """Print the memory usage information."""
    if verbose:
        print(
            f"Memory usage decreased from {start_mem:.2f} MB to {end_mem:.2f} MB "
            f"({100 * (start_mem - end_mem) / start_mem:.1f}% reduction)"
        )


def get_column_name_mappings():
    """Configures the column name mappings between different DIA search software and DIA-NN.

    This function defines a dictionary that maps column names from Spectronaut and OpenSWATH
    outputs to the corresponding DIA-NN column names.

    Returns:
        dict: A dictionary with software names as keys and mapping dictionaries as values.
    """
    mappings = {
        "Spectronaut": {
            "ProteinName": "Protein.Ids",
            "PeptideSequence": "Stripped.Sequence",
            # Add more Spectronaut-specific mappings here...
        },
        "OpenSWATH": {
            "protein_id": "Protein.Ids",
            "sequence": "Stripped.Sequence",
            # Add more OpenSWATH-specific mappings here...
        }
        # Add other software mappings if necessary
    }
    return mappings


def harmonize_columns(dataframe, software):
    """
    Harmonizes the column names of a DataFrame to match DIA-NN naming conventions.

    Depending on the specified software, it renames the columns in the DataFrame
    to align with the DIA-NN column naming convention. It is intended for use with
    DataFrame objects containing data from either Spectronaut or OpenSWATH DIA search
    algorithms. If the data is already from DIA-NN, no renaming is performed.

    Args:
        dataframe (pd.DataFrame): The DataFrame with the original column names.
        software (str): The software name from which the data originated.
                        Expected values are 'DIA-NN', 'Spectronaut', or 'OpenSWATH'.

    Returns:
        pd.DataFrame: A DataFrame with harmonized column names.

    Raises:
        ValueError: If the software name is not recognized or mappings are not defined.

    Examples:
        >>> # Assuming you have dataframes df_spectronaut and df_openswath from the respective software.
        >>> df_spectronaut_harmonized = harmonize_columns(df_spectronaut, 'Spectronaut')
        >>> df_openswath_harmonized = harmonize_columns(df_openswath, 'OpenSWATH')
        >>> # The dataframes now have columns named according to DIA-NN conventions.
    """
    # Retrieve the column name mappings using the configuration function
    column_name_mappings = get_column_name_mappings()

    if software == "DIA-NN":
        return dataframe  # If it's already DIA-NN, no changes are needed

    if software in column_name_mappings:
        # Get the specific mapping for the provided software
        column_mapping = column_name_mappings[software]

        # Identify columns present in both the DataFrame and the mapping, then rename them
        common_columns = set(dataframe.columns) & set(column_mapping.keys())
        dataframe = dataframe.rename(
            columns={col: column_mapping[col] for col in common_columns}
        )
    else:
        raise ValueError(
            f"Column name mappings for software '{software}' are not defined."
        )

    return dataframe


def assess_fdr(
    dataframe,
    peptide_fdr_column,
    protein_fdr_column,
    pep_fdr_threshold=0.01,
    pro_fdr_threshold=0.01,
):
    """
    Performs a two-dimensional False Discovery Rate (FDR) assessment on proteomic data.

    This function scrutinizes false discoveries at both peptide and protein levels by applying
    the specified FDR threshold. It filters the data to retain only those identifications
    that meet the FDR criteria at both levels, enhancing the precision and reliability of the findings.

    Args:
        dataframe (pd.DataFrame): DataFrame containing proteomic data.
        peptide_fdr_column (str): Column name in the DataFrame for peptide-level FDR values.
        protein_fdr_column (str): Column name in the DataFrame for protein-level FDR values.
        pep_fdr_threshold (float, optional): FDR threshold for filtering the peptides (default is 0.01).
        pro_fdr_threshold (float, optional): FDR threshold for filtering the proteins (default is 0.01).

    Returns:
        pd.DataFrame: A filtered DataFrame containing only the rows where both peptide-level
                    and protein-level FDR values are below the threshold.

    Examples:
        >>> proteomics_data = pd.DataFrame({
                'Peptide': ['pep1', 'pep2', 'pep3'],
                'Protein': ['prot1', 'prot2', 'prot3'],
                'Peptide_FDR': [0.005, 0.02, 0.01],
                'Protein_FDR': [0.01, 0.005, 0.02]
            })
        >>> filtered_data = assess_fdr(proteomics_data, 'Peptide_FDR', 'Protein_FDR')
        >>> print(filtered_data)
        # This will display rows from proteomics_data where both Peptide_FDR and Protein_FDR are below 0.01
    """

    # Filter the DataFrame based on the FDR thresholds
    filtered_df = dataframe[
        (dataframe[peptide_fdr_column] < pep_fdr_threshold)
        & (dataframe[protein_fdr_column] < pro_fdr_threshold)
    ]

    return filtered_df


def data_transform(df, convert_reverse=False) -> pd.DataFrame:
    """
    Perform index transformation on the input DataFrame.

    This function either combines multiple columns into a single index or splits an existing
    index into multiple columns, depending on the `convert_reverse` flag. When `convert_reverse`
    is False, it concatenates specified columns to create a new index. When True, it splits the
    index into separate columns.

    Parameters:
    ----------
    df : pd.DataFrame
        The input DataFrame to be transformed.
    convert_reverse : bool, optional
        If True, performs the reverse operation: splits the index into separate columns.
        Default is False, which concatenates specified columns to form a new index.

    Returns:
    -------
    pd.DataFrame
        The transformed DataFrame with either a new combined index or split columns
        from the index, based on the `convert_reverse` flag.

    Examples:
    --------
    >>> df = pd.DataFrame({
            'Genes': ['Gene1', 'Gene2'],
            'Stripped.Sequence': ['Seq1', 'Seq2'],
            'Modified.Sequence': ['Mod1', 'Mod2'],
            'Precursor.Id': ['Prec1', 'Prec2'],
            'file_name': ['File1', 'File2']
        })
    >>> transformed_df = data_transform(df)
    >>> print(transformed_df.index)
    # This will display a new index created by concatenating the specified columns.

    >>> reverse_transformed_df = data_transform(df, convert_reverse=True)
    >>> print(reverse_transformed_df.columns)
    # This will display the columns split from the original index.
    """
    if convert_reverse:
        df["index"] = df.index
        df = df.loc[:, ~df.columns.duplicated()]
        df[
            [
                "Protein.Group",
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
                "Protein.Group",
                "Stripped.Sequence",
                "Modified.Sequence",
                "Precursor.Id",
                "file_name",
            ],
            inplace=True,
        )
    else:
        df["index"] = (
            df["Protein.Group"]
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


def construct_index(df, convert_reverse=False) -> pd.DataFrame:
    """
    Perform index transformation on the input DataFrame.

    This function either combines multiple columns into a single index or splits an existing
    index into multiple columns, depending on the `convert_reverse` flag. When `convert_reverse`
    is False, it concatenates specified columns to create a new index. When True, it splits the
    index into separate columns.

    Parameters:
    ----------
    df : pd.DataFrame
        The input DataFrame to be transformed.
    convert_reverse : bool, optional
        If True, performs the reverse operation: splits the index into separate columns.
        Default is False, which concatenates specified columns to form a new index.

    Returns:
    -------
    pd.DataFrame
        The transformed DataFrame with either a new combined index or split columns
        from the index, based on the `convert_reverse` flag.

    Examples:
    --------
    >>> df = pd.DataFrame({
            'Genes': ['Gene1', 'Gene2'],
            'Stripped.Sequence': ['Seq1', 'Seq2'],
            'Modified.Sequence': ['Mod1', 'Mod2'],
            'Precursor.Id': ['Prec1', 'Prec2'],
            'file_name': ['File1', 'File2']
        })
    >>> transformed_df = construct_index(df)
    >>> print(transformed_df.index)
    # This will display a new index created by concatenating the specified columns.

    >>> reverse_transformed_df = construct_index(df, convert_reverse=True)
    >>> print(reverse_transformed_df.columns)
    # This will display the columns split from the original index.
    """
    if convert_reverse:
        df = df.assign(
            **{
                "Protein.Group": df.index.get_level_values(0),
                "Stripped.Sequence": df.index.get_level_values(1),
                "Modified.Sequence": df.index.get_level_values(2),
                "Precursor.Id": df.index.get_level_values(3),
                "file_name": df.index.get_level_values(4),
            }
        )
        df.set_index(
            [
                "Protein.Group",
                "Stripped.Sequence",
                "Modified.Sequence",
                "Precursor.Id",
                "file_name",
            ],
            inplace=True,
        )
    else:
        df["index"] = (
            df["Protein.Group"]
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


def coefficient_variation(data, log_transformed=False, dropna_thresh=0.5):
    """
    Calculate the coefficient of variation for each plasma protein in the dataset.

    This function computes the standard deviation divided by the mean (coefficient of variation)
    for each row in the data. It can optionally apply a log10 transformation to the data before
    the calculation. Rows with a significant number of NaN values (as determined by dropna_thresh)
    are dropped from the calculation.

    Args:
        data (pd.DataFrame): DataFrame containing the plasma protein data.
        log_transformed (bool): If True, applies a log2 transformation to the data. Defaults to False.
        dropna_thresh (float): Threshold for dropping rows with NaN values. Defaults to 0.5.

    Returns:
        pd.DataFrame: DataFrame containing the original data along with the calculated coefficient
                      of variation, mean abundance, and abundance rank for each protein.

    Examples:
        >>> data = pd.DataFrame({
                'Protein1': [1, 2, 3, np.nan],
                'Protein2': [4, 5, np.nan, 7]
            })
        >>> result = coefficient_variation(data)
        >>> print(result)
    """
    if log_transformed:
        data = np.log2(data)

    data = data.dropna(axis=0, thresh=len(data.columns) * dropna_thresh)
    mean_abundance = data.mean(axis=1)
    std_deviation = data.std(axis=1)
    data["Coefficient of Variation [%]"] = (std_deviation / mean_abundance) * 100
    data["Abundance"] = mean_abundance
    data = data.sort_values("Abundance", ascending=False)
    data["Abundance_Rank"] = range(1, len(data) + 1)

    return data
