#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : MaxLFQ.py
@Time     : 2022/11/11 12:55:03
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
'''
# here put the import lib

import numpy as np
import pandas as pd
import numba

from scipy.linalg import lstsq


def maxLFQ(X):
    """
    Perform the maxLFQ algorithm on a given set of protein/peptide abundance data.

    This function applies the maxLFQ (maximum label-free quantification) algorithm,
    which is used for estimating the abundance of proteins/peptides in mass spectrometry-based
    proteomics data. It includes handling of missing values and grouping of samples for comparison.

    Args:
        X (pd.DataFrame): A DataFrame containing the protein/peptide abundance data. Rows represent
                            proteins/peptides, and columns represent samples.

    Returns:
        dict: A dictionary with two keys: 'estimate' containing the estimated abundances and
                'annotation' containing group annotations.

    Example:
        >>> data = pd.DataFrame(np.random.rand(10, 5), columns=['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'])
        >>> result = maxLFQ(data)
        >>> print(result)

    Note:
        This is an optimized version of the algorithm, designed for large-scale computations.
    """
    try:
        # Check if DataFrame is empty or all values are NaN
        if X.empty:
            raise ValueError("Input DataFrame is empty.")
        if X.isna().all().all():
            raise ValueError("All values in the input DataFrame are NaN.")

        # Handle single row case
        if len(X) == 1:
            return pd.DataFrame({'estimate': X.iloc[0].tolist(), 'annotation': [""]})
        
        N = X.shape[1]

        def spread(i, g, cc, X):
            # Spread function for assigning group labels
            g[i] = cc
            for r in range(X.shape[0]):
                if not pd.isna(X.iloc[r, i]):
                    for k in range(X.shape[1]):
                        if not pd.isna(X.iloc[r, k]) and pd.isna(g[k]):
                            spread(k, g, cc, X)

        def maxLFQ_do(X_sub):
            # maxLFQ algorithm implementation
            N = X_sub.shape[1]
            AtA = np.zeros((N, N))
            Atb = np.zeros(N)

            for i in range(N - 1):
                for j in range(i + 1, N):
                    r_i_j = np.nanmedian(X_sub.iloc[:, j] - X_sub.iloc[:, i])
                    if not np.isnan(r_i_j):
                        AtA[i, j] = AtA[j, i] = -1
                        AtA[i, i] += 1
                        AtA[j, j] += 1
                        Atb[i] -= r_i_j
                        Atb[j] += r_i_j

            A = np.vstack([np.hstack([2 * AtA, np.ones(N).reshape(-1, 1)]), 
                           np.hstack([np.ones(N), np.zeros(1)])])
            b = np.append(2 * Atb, np.nanmean(X_sub.to_numpy()) * N)

            res = lstsq(A, b, lapack_driver='gelsy')[0]
            return res[:N]

        # Applying the maxLFQ algorithm
        cc = 0
        g = [np.nan] * N
        for i in range(N):
            if pd.isna(g[i]):
                cc += 1
                spread(i, g, cc, X)

        w = [np.nan] * N
        for i in range(1, cc + 1):
            ind = [idx for idx, val in enumerate(g) if val == i]
            if len(ind) == 1:
                w[ind[0]] = np.nanmedian(X.iloc[:, ind[0]])
            else:
                X_sub = X.iloc[:, ind]
                results = maxLFQ_do(X_sub)
                for index, res in zip(ind, results):
                    w[index] = res

        # Check for all NaN results
        if np.isnan(w).all():
            # return {"estimate": w, "annotation": "NA"}
            raise ValueError("maxLFQ calculation resulted in all NaN values.")

        # Preparing final results
        annotation = ";".join(map(str, g)) if np.any(np.isnan(w)) else ""
        return {"estimate": w, "annotation": annotation}

    except ValueError as ve:
        return {"error": str(ve)}
    except Exception as e:
        return {"error": f"An unexpected error occurred: {str(e)}"}


def generate_peptide_data(n_peptides, n_samples, missing_data_rate=0.1):
    """
    Generates proteomics data with missing values.

    This function creates a DataFrame of proteomics data (e.g., peptide or protein abundances)
    for a given number of peptides and samples. It introduces missing values to simulate a
    typical proteomics dataset.

    Args:
        n_peptides (int): The number of peptides (or proteins), corresponding to the number of rows.
        n_samples (int): The number of samples, corresponding to the number of columns.
        missing_data_rate (float, optional): The proportion of missing values in each sample. Defaults to 0.1.

    Returns:
        pd.DataFrame: A DataFrame containing the generated proteomics data. Each row represents a peptide (or protein),
                      and each column represents a sample. Missing data are represented by NaN values.

    Example:
        >>> peptide_data = generate_peptide_data(100, 10, 0.1)
        >>> print(peptide_data.head())
    """
    np.random.seed(0)
    data = np.random.lognormal(mean=2, sigma=0.5, size=(n_peptides, n_samples))

    # Introduce missing values, ensuring approximately 'missing_data_rate' missingness per sample
    for col in range(n_samples):
        missing_indices = np.random.choice(n_peptides, size=int(n_peptides * missing_data_rate), replace=False)
        data[missing_indices, col] = np.nan

    return pd.DataFrame(data, columns=[f'Sample_{i+1}' for i in range(n_samples)])


def create_synthetic_test_data(num_proteins=100, num_samples=10, max_peptides_per_protein=10):
    """
    Creates a synthetic DataFrame mimicking DIA proteomics data.

    This function generates synthetic test data for proteomics analysis. It simulates
    multiple proteins, each with a random number of peptides, and assigns random
    intensity values to these peptides across different samples. Missing values are
    introduced to simulate a realistic DIA dataset.

    Args:
        num_proteins (int): The number of proteins to simulate. Defaults to 100.
        num_samples (int): The number of samples for each peptide. Defaults to 10.
        max_peptides_per_protein (int): The maximum number of peptides per protein. Defaults to 10.

    Returns:
        pd.DataFrame: A DataFrame representing synthetic proteomics data. Rows correspond to peptides,
                      columns to samples, and values to peptide intensities.

    Example:
        >>> synthetic_data = create_synthetic_test_data(num_proteins=100, num_samples=10, max_peptides_per_protein=10)
        >>> print(synthetic_data.head())
    """
    data = []
    peptide_info = []  # Stores information about peptides and proteins

    for protein_id in range(num_proteins):
        num_peptides = np.random.randint(1, max_peptides_per_protein + 1)
        for peptide_id in range(num_peptides):
            mean_intensity = np.random.uniform(1, 3)  # Random average intensity
            sigma_intensity = np.random.uniform(0.3, 1.0)  # Random standard deviation
            peptide_intensity = np.random.lognormal(mean=mean_intensity, sigma=sigma_intensity, size=num_samples)

            # Introduce missing values more likely in low-intensity peptides
            missing_probability = np.clip(1 / peptide_intensity, 0, 1)
            missing_indices = np.random.uniform(size=num_samples) < missing_probability
            peptide_intensity[missing_indices] = np.nan

            data.append(peptide_intensity)
            peptide_info.append(f"Protein_{protein_id}_Peptide_{peptide_id}")

    # Introduce sample-specific perturbations
    sample_specific_effects = np.random.normal(loc=0, scale=0.2, size=num_samples)
    data = np.array(data) + sample_specific_effects

    # Create DataFrame
    df = pd.DataFrame(data.transpose(), columns=peptide_info)
    df.index.name = 'Sample'
    return df.T


def maxLFQ_optimal(X):
    """
    Perform the maxLFQ algorithm on a given set of protein/peptide abundance data.

    This function applies the maxLFQ (maximum label-free quantification) algorithm,
    which is used for estimating the abundance of proteins/peptides in mass spectrometry-based
    proteomics data. It includes handling of missing values and grouping of samples for comparison.

    Args:
        X (np.ndarray): A 2D array containing the protein/peptide abundance data. Rows represent
                        proteins/peptides, and columns represent samples.

    Returns:
        dict: A dictionary with two keys: 'estimate' containing the estimated abundances and
                'annotation' containing group annotations.

    Example:
        >>> data = np.random.rand(10, 5)
        >>> result = maxLFQ(data)
        >>> print(result)

    Note:
        This is an optimized version of the algorithm, designed for large-scale computations.
    """
    try:
        # Check if array is empty or all values are NaN
        if X.size == 0:
            raise ValueError("Input array is empty.")
        if np.isnan(X).all():
            raise ValueError("All values in the input array are NaN.")

        # Handle single row case
        if X.shape[0] == 1:
            return {'estimate': X[0], 'annotation': [""]}

        N = X.shape[1]

        def spread(i, g, cc):
            # Spread function for assigning group labels
            g[i] = cc
            for r in range(X.shape[0]):
                if not np.isnan(X[r, i]):
                    for k in range(X.shape[1]):
                        if not np.isnan(X[r, k]) and np.isnan(g[k]):
                            spread(k, g, cc)

        def maxLFQ_do(X_sub):
            # maxLFQ algorithm implementation
            N = X_sub.shape[1]
            AtA = np.zeros((N, N))
            Atb = np.zeros(N)

            for i in range(N - 1):
                for j in range(i + 1, N):
                    r_i_j = np.nanmedian(X_sub[:, j] - X_sub[:, i])
                    if not np.isnan(r_i_j):
                        AtA[i, j] = AtA[j, i] = -1
                        AtA[i, i] += 1
                        AtA[j, j] += 1
                        Atb[i] -= r_i_j
                        Atb[j] += r_i_j

            A = np.vstack([np.hstack([2 * AtA, np.ones((N, 1))]), 
                           np.hstack([np.ones(N), 0])])
            b = np.append(2 * Atb, np.nanmean(X_sub) * N)

            res = lstsq(A, b, lapack_driver='gelsy')[0]
            return res[:N]

        # Applying the maxLFQ algorithm
        cc = 0
        g = np.full(N, np.nan)
        for i in range(N):
            if np.isnan(g[i]):
                cc += 1
                spread(i, g, cc)

        w = np.full(N, np.nan)
        for i in range(1, cc + 1):
            ind = np.where(g == i)[0]
            if len(ind) == 1:
                w[ind[0]] = np.nanmedian(X[:, ind[0]])
            else:
                X_sub = X[:, ind]
                results = maxLFQ_do(X_sub)
                w[ind] = results

        # Check for all NaN results
        if np.isnan(w).all():
            raise ValueError("maxLFQ calculation resulted in all NaN values.")

        # Preparing final results
        annotation = ";".join(map(str, g)) if np.any(np.isnan(w)) else ""
        return {"estimate": w, "annotation": annotation}

    except ValueError as ve:
        return {"error": str(ve)}
    except Exception as e:
        return {"error": f"An unexpected error occurred: {str(e)}"}


@numba.jit
def spread(i, g, cc, X):
    """
    Recursively assigns group labels to elements in a matrix.

    Args:
        i (int): The current column index in the matrix.
        g (np.ndarray): An array to store group labels.
        cc (int): The current group label.
        X (np.ndarray): The input data matrix.

    This function is part of the maxLFQ algorithm and is optimized with Numba.
    It assigns group labels to elements in X, based on their non-NaN status and connectivity.
    """
    g[i] = cc
    for r in range(X.shape[0]):
        if not np.isnan(X[r, i]):
            for k in range(X.shape[1]):
                if not np.isnan(X[r, k]) and np.isnan(g[k]):
                    spread(k, g, cc, X)

@numba.jit
def maxLFQ_do(X_sub):
    """
    Performs the maxLFQ algorithm on a subset of the data.

    Args:
        X_sub (np.ndarray): A subset of the original data matrix.

    Returns:
        np.ndarray: The results of the maxLFQ algorithm on X_sub.

    This function implements the maxLFQ algorithm, optimized with Numba. It calculates the relative
    quantities of peptides/proteins from the input LC-MS/MS data.
    """
    N = X_sub.shape[1]
    AtA = np.zeros((N, N))
    Atb = np.zeros(N)

    for i in range(N - 1):
        for j in range(i + 1, N):
            r_i_j = np.nanmedian(X_sub[:, j] - X_sub[:, i])
            if not np.isnan(r_i_j):
                AtA[i, j] = AtA[j, i] = -1
                AtA[i, i] += 1
                AtA[j, j] += 1
                Atb[i] -= r_i_j
                Atb[j] += r_i_j

    A = np.vstack([np.hstack([2 * AtA, np.ones((N, 1))]),
                    np.hstack([np.ones(N), 0])])
    b = np.append(2 * Atb, np.nanmean(X_sub) * N)

    res = lstsq(A, b, lapack_driver='gelsy')[0]

    return res[:N]

def fast_maxLFQ(X):
    """
    Applies the maxLFQ algorithm to the given data matrix.

    Args:
        X (np.ndarray): A data matrix where rows represent samples and columns represent features.

    Returns:
        dict: A dictionary with two keys: 'estimate' containing the quantification result, and
              'annotation' containing information about the calculation or errors.

    This function checks for empty or NaN-only data, handles single-row data, assigns group labels,
    and applies the maxLFQ algorithm to each group of features.

    Usage:
        >>> data = np.array([[np.nan, 1, 2], [np.nan, np.nan, 2], [3, 4, np.nan]])
        >>> result = fast_maxLFQ(data)
        >>> print(result)
        {'estimate': [resulting values], 'annotation': [annotation string]}
    """
    if X.size == 0:
        return {"estimate": np.nan, "annotation": "Empty array"}
    if np.isnan(X).all():
        return {"estimate": np.nan, "annotation": "All NaN values"}

    # 处理单行数据
    if X.shape[0] == 1:
        return {"estimate": X[0], "annotation": [""]}

    N = X.shape[1]

    # Accelerated function calls using numba
    cc = 0
    g = np.full(N, np.nan)
    for i in range(N):
        if np.isnan(g[i]):
            cc += 1
            spread(i, g, cc, X)

    w = np.full(N, np.nan)
    for i in range(1, cc + 1):
        ind = np.where(g == i)[0]
        if len(ind) == 1:
            w[ind[0]] = np.nanmedian(X[:, ind[0]])
        else:
            X_sub = X[:, ind]
            results = maxLFQ_do(X_sub)
            w[ind] = results

    # Check for all NaN results
    if np.isnan(w).all():
        return {"estimate": w, "annotation": "NA"}
    
    # Preparing final results
    annotation = ";".join(map(str, g)) if np.any(np.isnan(w)) else ""
    return {"estimate": w, "annotation": annotation}


class MaxLFQOptimizer:
    """
    A class for optimizing the maxLFQ algorithm using Numba for faster computation.

    This class contains methods for performing the maxLFQ quantification algorithm on mass spectrometry data.
    It is optimized with Numba to improve performance.

    Methods:
        spread(i, g, cc, X): Assigns group labels to elements in a matrix.
        maxLFQ_do(X_sub): Performs the maxLFQ algorithm on a subset of the data.
        maxLFQ_optimal_v2(X): Applies the maxLFQ algorithm to the given data matrix.

    Usage:
        >>> optimizer = MaxLFQOptimizer()
        >>> data = np.array([[np.nan, 1, 2], [np.nan, np.nan, 2], [3, 4, np.nan]])
        >>> result = optimizer.maxLFQ_optimal_v2(data)
        >>> print(result)

    Note:
        This is an optimized version of the algorithm, designed for large-scale computations.
    """

    def __init__(self):
        """ Initializes the MaxLFQOptimizer class. """
        pass

    @staticmethod
    @numba.jit
    def _spread(i, g, cc, X):
        """
        Recursively assigns group labels to elements in a matrix.

        Args:
            i (int): The current column index in the matrix.
            g (np.ndarray): An array to store group labels.
            cc (int): The current group label.
            X (np.ndarray): The input data matrix.

        This method is part of the maxLFQ algorithm and is optimized with Numba.
        It assigns group labels to elements in X, based on their non-NaN status and connectivity.
        """
        g[i] = cc
        for r in range(X.shape[0]):
            if not np.isnan(X[r, i]):
                for k in range(X.shape[1]):
                    if not np.isnan(X[r, k]) and np.isnan(g[k]):
                        MaxLFQOptimizer._spread(k, g, cc, X)

    @staticmethod
    @numba.jit
    def _maxLFQ_do(X_sub):
        """
        Performs the maxLFQ algorithm on a subset of the data.

        Args:
            X_sub (np.ndarray): A subset of the original data matrix.

        Returns:
            np.ndarray: The results of the maxLFQ algorithm on X_sub.

        This method implements the maxLFQ algorithm, optimized with Numba. It calculates the relative
        quantities of peptides/proteins from the input LC-MS/MS data.
        """
        N = X_sub.shape[1]
        AtA = np.zeros((N, N))
        Atb = np.zeros(N)

        for i in range(N - 1):
            for j in range(i + 1, N):
                r_i_j = np.nanmedian(X_sub[:, j] - X_sub[:, i])
                if not np.isnan(r_i_j):
                    AtA[i, j] = AtA[j, i] = -1
                    AtA[i, i] += 1
                    AtA[j, j] += 1
                    Atb[i] -= r_i_j
                    Atb[j] += r_i_j

        A = np.vstack([np.hstack([2 * AtA, np.ones((N, 1))]),
                        np.hstack([np.ones(N), 0])])
        b = np.append(2 * Atb, np.nanmean(X_sub) * N)

        res = lstsq(A, b, lapack_driver='gelsy')[0]

        return res[:N]

    def maxLFQ_fast(self, X):
        """
        Applies the maxLFQ algorithm to the given data matrix.

        Args:
            X (np.ndarray): A data matrix where rows represent samples and columns represent features.

        Returns:
            dict: A dictionary with two keys: 'estimate' containing the quantification result, and
                  'annotation' containing information about the calculation or errors.

        This method checks for empty or NaN-only data, handles single-row data, assigns group labels,
        and applies the maxLFQ algorithm to each group of features.
        """
        if X.size == 0:
            return {"estimate": np.nan, "annotation": "Empty array"}
        if np.isnan(X).all():
            return {"estimate": np.nan, "annotation": "All NaN values"}

        # Handle single row case
        if X.shape[0] == 1:
            return {"estimate": X[0], "annotation": [""]}

        N = X.shape[1]

        # 使用 numba 加速的函数调用
        cc = 0
        g = np.full(N, np.nan)
        for i in range(N):
            if np.isnan(g[i]):
                cc += 1
                self._spread(i, g, cc, X)

        w = np.full(N, np.nan)
        for i in range(1, cc + 1):
            ind = np.where(g == i)[0]
            if len(ind) == 1:
                w[ind[0]] = np.nanmedian(X[:, ind[0]])
            else:
                X_sub = X[:, ind]
                results = self._maxLFQ_do(X_sub)
                w[ind] = results

        # Check for all NaN results
        if np.isnan(w).all():
            return {"estimate": w, "annotation": "NA"}
        
        # Preparing final results
        annotation = ";".join(map(str, g)) if np.any(np.isnan(w)) else ""
        return {"estimate": w, "annotation": annotation}


if __name__ == '__main__':
    # Example usage
    # data = pd.DataFrame(np.random.rand(10, 5), columns=['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'])
    # result = maxLFQ(data)
    # print(result)
    peptide_intensities1 = generate_peptide_data(10, 1000, 0.22)
    estinmate_intensities = maxLFQ(peptide_intensities1)
    print(estinmate_intensities)
