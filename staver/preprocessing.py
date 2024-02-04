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
import pymc3 as pm
import arviz as az
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import os

from joblib import Parallel, delayed


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


class PeptideBayesianEstimator:
    """
    A class for estimating peptide abundance using Bayesian statistics.

    Attributes:
        original_data (np.array): The original peptide data, including NaN values.
        data (np.array): The cleaned peptide data, NaN values removed.

    Examples:
        >>> data = np.array([1, 2, 3, 4, 5, np.nan, 7, 8, 9, 10])
        >>> peptide_data = PeptideBayesianEstimator(data)
        >>> trace = peptide_data.estimate_parameters(n_cores=4)
        >>> estimates = peptide_data.summarize_posterior(trace)
        >>> statistics = peptide_data.compute_observed_statistics()
        >>> peptide_data.plot_density(trace)
        >>> peptide_data.model_diagnostics(trace)
    """

    def __init__(self, data):
        """
        Initialize the PeptideData object with the provided data.

        Args:
        data (np.array): The peptide data.
        """
        self.original_data = data
        self.data = data[~np.isnan(data)]  # Drop rows with NaN values
        # self.identification_frequency = self.compute_identification_frequency(data)

    def estimate_parameters(
        self, use_variational_inference=True, n_samples=500, n_cores=4
    ):
        """
        Estimate parameters using a Bayesian model with options for variational inference or MCMC.
        Args:
            use_variational_inference (bool): Flag to use variational inference instead of MCMC.
            n_samples (int): Number of posterior samples to draw.
            n_cores (int): Number of cores for parallel processing.
        Returns:
            MultiTrace or InferenceData: A trace object containing the posterior distributions.
        """
        try:
            with pm.Model() as self.model:
                # Define the model (with simplified parameters if applicable)

                # Define priors
                mu = pm.Normal("mu", mu=np.mean(self.data), sd=np.std(self.data))
                sigma = pm.HalfNormal("sigma", sd=np.std(self.data))

                # Define likelihood
                likelihood = pm.Normal("y", mu=mu, sd=sigma, observed=self.data)

                if use_variational_inference:
                    # Using variational inference (ADVI) for faster estimation
                    logging.info("Using variational inference for estimation.")
                    approx = pm.fit(n=20000, method="advi")
                    trace = approx.sample(draws=n_samples)
                else:
                    # Using MCMC with optimized initial values and efficient sampler
                    logging.info("Using MCMC for estimation.")
                    start = pm.find_MAP()
                    step = pm.NUTS()
                    trace = pm.sample(
                        draws=n_samples, tune=500, step=step, cores=n_cores, start=start
                    )

            return trace
        except Exception as e:
            logging.error("Error in parameter estimation: %s", e)
            raise

    # def run_inference(self, draws=3000, tune=1000):
    #     """
    #     Run Bayesian inference using MCMC.

    #     Parameters:
    #     draws (int): Number of samples to draw from the posterior distribution.
    #     tune (int): Number of steps to tune the sampler.

    #     Returns:
    #     trace (MultiTrace): The MCMC trace from the model.
    #     """
    #     model = self.estimate_parameters()
    #     with model:
    #         # Run MCMC
    #         trace = pm.sample(draws, tune=tune, return_inferencedata=False)

    #     return trace

    def summarize_posterior(self, trace):
        """
        Summarize the posterior distribution of the parameters.

        Args:
            trace (MultiTrace): The MCMC trace from the model.

        Returns:
            DataFrame: Summary statistics of the posterior distribution.
        """
        return az.summary(trace)

    def compute_identification_frequency(self, data):
        """
        Compute the identification frequency for each peptide.
        Args:
            data (np.array): The peptide data.
        Returns:
            np.array: The identification frequencies.
        """
        if data.ndim == 1:
            # 处理一维数组的情况
            return np.count_nonzero(~np.isnan(data)) / data.size
        elif data.ndim == 2:
            # 处理二维数组的情况
            return np.count_nonzero(~np.isnan(data), axis=1) / data.shape[1]
        else:
            raise ValueError("Data must be either 1D or 2D array")

    def compute_posterior_statistics(self, trace, credible_interval=0.95):
        """
        Compute the posterior summary statistics for each parameter in the trace.

        Args:
            trace (MultiTrace): The MCMC trace from the model.
            credible_interval (float): The credible interval for the summary statistics.

        Returns:
            DataFrame: A dataframe containing the posterior summary statistics.
        """
        # 使用 ArviZ 计算后验分布的总结统计数据
        summary = az.summary(trace, hdi_prob=credible_interval)

        # 计算每个肽段的鉴定频率
        identification_frequency = self.compute_identification_frequency(
            self.original_data
        )
        summary["identification_frequency"] = identification_frequency

        # 计算变异系数 (CV)
        for var in trace.varnames:
            if not var.endswith("__"):  # 忽略 PyMC3 内部变量
                mean = summary.loc[var, "mean"]
                sd = summary.loc[var, "sd"]
                summary.loc[var, "cv"] = sd / mean if mean != 0 else np.nan

        return summary

    def compute_observed_statistics(self):
        """
        Compute the identification frequency, mean, median, and coefficient of variation (CV).

        Returns:
            dict: A dictionary containing the computed statistics.
        """
        identification_frequency = self.compute_identification_frequency(
            self.original_data
        )
        mean_abundance = np.mean(self.data)
        median_abundance = np.median(self.data)
        standard_deviation = np.std(self.data)
        cv = stats.variation(self.data)

        return {
            "Identification Frequency": identification_frequency,
            "Mean Abundance": mean_abundance,
            "Median Abundance": median_abundance,
            "Standard Deviation": standard_deviation,
            "Coefficient of Variation": cv,
        }

    def plot_density(self, trace):
        """
        Plot density graphs for observed data and posterior predictions.

        Args:
            trace (MultiTrace): The MCMC trace from the model.
        """
        # 为观测值绘制密度图
        sns.kdeplot(self.data, color="blue", label="Observed Data")

        # 从后验分布生成预测数据并绘制密度图
        post_pred = pm.sample_posterior_predictive(trace, model=self.model, samples=500)
        post_pred_data = post_pred["y"].flatten()
        sns.kdeplot(post_pred_data, color="red", label="Posterior Predictions")

        plt.xlabel("Peptide Abundance")
        plt.ylabel("Density")
        plt.title("Density Plot of Peptide Data")
        plt.legend()
        plt.show()

    def model_diagnostics(self, trace):
        """
        Perform model diagnostics, including MCMC convergence checks and posterior predictive checks.

        Args:
            trace (MultiTrace): The MCMC trace from the model.
        """
        # MCMC收敛诊断
        # Gelman-Rubin 统计量
        gelman_rubin = az.rhat(trace)
        print("Gelman-Rubin Diagnostic:\n", gelman_rubin)

        # Posterior predictive checking
        with self.model:
            observed_var_name = "y"
            ppc = pm.sample_posterior_predictive(trace, var_names=[observed_var_name])

        # Plotting a posteriori predictions against actual observations
        az.plot_ppc(az.from_pymc3(posterior_predictive=ppc, model=self.model))
        plt.show()


def extract_feature_stats(df):
    """
    Extract feature statistics for mean abundance, standard deviation, and identification frequency.

    Args:
        df (pd.DataFrame): The DataFrame containing the summary statistics from Bayesian analysis.

    Returns:
        dict: A dictionary with extracted statistics.
    """
    # 提取统计值
    mean_abundance = df.loc["mu", "mean"]
    standard_deviation = df.loc["sigma", "mean"]
    identification_frequency = df.loc["mu", "identification_frequency"]

    return {
        "mean_abundance": mean_abundance,
        "standard_deviation": standard_deviation,
        "identification_frequency": identification_frequency,
    }


def process_peptide(row, n_iterations, n_cores):
    # Function to process each peptide
    estimator = PeptideBayesianEstimator(row)
    trace = estimator.estimate_parameters(n_iterations, n_cores=n_cores)
    return estimator.compute_posterior_statistics(trace)


def estimate_peptide_features(dataframe, n_iterations=1000, n_cores=10):
    """
    Estimate features for each peptide in the dataframe using Bayesian statistics.

    Args:
        dataframe (pd.DataFrame): DataFrame containing peptide data.
        n_iterations (int): Number of iterations for MCMC sampling.
        n_cores (int): Number of cores for parallel sampling.

    Returns:
        pd.DataFrame: A DataFrame with estimated features for each peptide.

    Examples:
        >>> df = pd.DataFrame(np.random.rand(10, 5))
        >>> estimated_features = estimate_peptide_features(df, n_iterations=500, n_cores=2)
        >>> print(estimated_features.head())
    """

    # Parallel processing with joblib
    results = Parallel(n_jobs=n_cores)(
        delayed(process_peptide)(row.to_numpy(), n_iterations, 1)
        for _, row in dataframe.iterrows()
    )

    # Post-process results to extract features and handle errors
    results = [
        extract_feature_stats(res) if res is not None else None for res in results
    ]

    # Convert results to DataFrame
    results_df = pd.DataFrame(results, index=dataframe.index)
    return results_df
