import pandas as pd
import numpy as np

from preprocessing import estimate_peptide_features


def calculate_distribution_metrics(dataframe):
    """
    Highly optimized function to calculate skewness and kurtosis for peptide abundance distribution using pure NumPy.

    Parameters:
    ----------
    dataframe : pd.DataFrame
        A DataFrame where each column represents a sample and each row represents peptide abundance.

    Returns:
    -------
    pd.DataFrame
        A DataFrame with calculated skewness and kurtosis for each peptide.
    """
    data = dataframe.to_numpy()

    # Mean and standard deviation
    mean = np.nanmean(data, axis=1, keepdims=True)
    std_dev = np.nanstd(data, axis=1, keepdims=True)

    # Skewness and kurtosis calculations using vectorized operations
    skewness = np.nanmean(((data - mean) / std_dev) ** 3, axis=1)
    kurtosis_vals = np.nanmean(((data - mean) / std_dev) ** 4, axis=1) - 3

    # Creating the result DataFrame
    distribution_metrics = pd.DataFrame(
        {"skewness": abs(skewness), "kurtosis": abs(kurtosis_vals)},
        index=dataframe.index,
    )

    return distribution_metrics


def remove_low_outliers(df, fraction=5):
    """
    Replace elements in the DataFrame that are below a fraction of the row mean with NaN using NumPy vectorization.

    Args:
    df (DataFrame): The input pandas DataFrame.
    fraction (int): The fraction of the row mean that elements must be greater than or equal to avoid replacement.

    Returns:
    DataFrame: A DataFrame with elements below the row mean fraction replaced with NaN.
    """
    data = df.values  # Convert DataFrame to NumPy array for faster computation
    row_means = np.nanmean(data, axis=1, keepdims=True)
    mask = data >= (row_means / fraction)
    return pd.DataFrame(
        np.where(mask, data, np.nan), index=df.index, columns=df.columns
    )


def intensity_stats(sampled_data):
    mean_intensity = sampled_data.mean(axis=1, skipna=True)
    median_intensity = sampled_data.median(axis=1, skipna=True)
    std_intensity = sampled_data.std(axis=1, skipna=True)
    cv_intensity = std_intensity / mean_intensity
    freq_peptide = sampled_data.count(axis=1) / sampled_data.shape[1]
    count_peptide = sampled_data.count(axis=1)

    stats = pd.DataFrame(
        {
            "mean_intensity": mean_intensity,
            "median_intensity": median_intensity,
            "std_intensity": std_intensity,
            "cv_intensity": cv_intensity,
            "frequency": freq_peptide,
            "count": count_peptide,
        }
    )
    return stats


def calculate_row_mad(df):
    """
    Calculate the Mean Absolute Deviation (MAD) for each row of the DataFrame using NumPy vectorization.

    Args:
    df (DataFrame): The input pandas DataFrame.

    Returns:
    Series: A pandas Series containing the MAD for each row.
    """
    cleaned_data = remove_low_outliers(df).values
    row_means = np.nanmean(cleaned_data, axis=1, keepdims=True)
    mad = np.nanmean(np.abs(cleaned_data - row_means), axis=1)
    return pd.Series(mad, index=df.index)


def retention_time_deviation(sampled_data):
    mean_RT = sampled_data.mean(axis=1, skipna=True)
    median_RT = sampled_data.median(axis=1, skipna=True)
    std_RT = sampled_data.std(axis=1, skipna=True)
    mad_RT = calculate_row_mad(sampled_data)
    shifted_RT = mad_RT / mean_RT

    rtd = pd.DataFrame(
        {
            "mean_RT": mean_RT,
            "median_RT": median_RT,
            "std_RT": std_RT,
            "mad_RT": mad_RT,
            "RTD": shifted_RT,
        }
    )
    return rtd


class PeptideSampling:
    """
    A class to perform bootstrap sampling on peptide intensity data.

    Attributes:
        data (pd.DataFrame): A dataframe with peptides as rows, samples as columns, and intensities as values.
        n_samples (int): The number of samples in the dataset.
        seed (int): Random seed for reproducibility.
        sampling_ratio (float): The ratio of samples to include in each bootstrap sample.
        n_bootstraps (int): The number of bootstrap samples to generate.

    Methods:
        set_seed(seed): Sets the random seed.
        set_sampling_ratio(ratio): Sets the sampling ratio.
        auto_adjust_bootstraps(): Automatically adjusts the number of bootstrap iterations.
        bootstrap(): Performs the bootstrap sampling.
        get_statistics(): Returns the bootstrap statistics as a dataframe.
    """

    def __init__(self, data, seed=42, sampling_ratio=0.8):
        """
        Initializes the PeptideSampling class with data, a random seed, and a sampling ratio.

        Args:
            data (pd.DataFrame): A dataframe with peptides as rows, samples as columns, and intensities as values.
            seed (int, optional): Random seed for reproducibility. Defaults to 42.
            sampling_ratio (float, optional): The ratio of samples to include in each bootstrap sample. Defaults to 0.5.
        """
        self.data = data
        self.n_samples = data.shape[1]
        self.seed = seed
        self.sampling_ratio = sampling_ratio
        self.n_bootstraps = self.auto_adjust_bootstraps()
        np.random.seed(self.seed)  # Set the random seed once for the entire process

    def set_sampling_ratio(self, ratio):
        """Sets the sampling ratio."""
        self.sampling_ratio = ratio

    def auto_adjust_bootstraps(self):
        """Automatically adjusts the number of bootstrap iterations based on the sample size."""
        return max(1000, 10 * self.n_samples)

    def bootstrap(self, calculate_intensity=True, calculate_RT=False):
        """
        Performs bootstrap sampling and records statistics for each peptide based on specified calculations.

        Args:
            calculate_intensity (bool, optional): Whether to calculate intensity statistics. Defaults to True.
            calculate_RT (bool, optional): Whether to calculate RT statistics. Defaults to False.
        """
        if not calculate_intensity and not calculate_RT:
            raise ValueError(
                "No statistical calculation selected. Choose at least one: intensity or RT."
            )

        bootstrap_stats = []

        for _ in range(self.n_bootstraps):
            sampled_cols = np.random.choice(
                self.data.columns,
                int(self.n_samples * self.sampling_ratio),
                replace=True,
            )
            sampled_data = self.data[sampled_cols]

            if calculate_intensity:
                # Calculating intensity statistics for each peptide
                stats = intensity_stats(sampled_data)
                # Calculating distribution metrics
                distribution_metrics = calculate_distribution_metrics(sampled_data)

            if calculate_RT:
                # Calculating RT statistics for each peptide
                stats = retention_time_deviation(sampled_data)

            bootstrap_stats.append(stats)

        self.bootstrap_stats = bootstrap_stats

    def get_statistics(self):
        """Returns the aggregated bootstrap statistics as a dataframe."""
        all_stats = pd.concat(self.bootstrap_stats)
        return all_stats, all_stats.groupby(level=0).mean()


def calculate_rt_consistency(df):
    tmp_wide_RT = df.pivot(index="index", columns="File.Name", values="Predicted.iRT")
    bootstrap = PeptideSampling(tmp_wide_RT, seed=42, sampling_ratio=0.8)
    bootstrap.bootstrap()
    _, RT_stats = bootstrap.get_statistics()
    return RT_stats["RTD"]


def extract_peptide_features(df):
    """
    Extracts peptide features from a DataFrame of peptide abundances.

    Args:
    dataframe (DataFrame): The input pandas DataFrame.

    Returns:
    DataFrame: A DataFrame containing the extracted peptide features.
    """
    bayes_estimates = estimate_peptide_features(df)
    bayes_estimates["CV"] = (
        bayes_estimates["standard_deviation"] / bayes_estimates["mean_abundance"]
    )
    bayes_estimates.rename(
        columns={
            "mean_abundance": "Abundance",
            "identification_frequency": "Frequency",
        },
        inplace=True,
    )
    distribution_metrics = calculate_distribution_metrics(df)
    rtd = calculate_rt_consistency(df)
    return pd.concat([bayes_estimates, distribution_metrics, rtd], axis=1)


def calculate_distribution_adjustment_factor(df):
    """
    Calculates the distribution adjustment factor for peptides.

    Args:
        df (pd.DataFrame): DataFrame containing 'skewness' and 'kurtosis' columns for peptides.

    Returns:
        pd.DataFrame: DataFrame with an additional column 'D_factor' representing the distribution adjustment factor for each peptide.
    """
    # Extract skewness and kurtosis
    skewness = df["skewness"]
    kurtosis = df["kurtosis"]

    # Calculate maximum absolute skewness and deviation of kurtosis from 3
    S_max = np.max(np.abs(skewness))
    K_max = np.max(np.abs(kurtosis - 3))

    # Calculate W_D
    skew_rsd = np.std(np.abs(skewness)) / np.mean(np.abs(skewness))
    kurt_rsd = np.std(np.abs(kurtosis)) / np.mean(np.abs(kurtosis))
    W_D = 1 / ((1 + (skew_rsd + kurt_rsd)) ** 3)

    print(f"S_max: {S_max}, K_max: {K_max}")
    print(f"skew_cv: {skew_rsd}; kurt_cv: {kurt_rsd}; W_D: {W_D}")
    # Calculate distribution adjustment factor D
    df["D_factor"] = (
        1 - ((np.abs(skewness) / S_max) + (np.abs(kurtosis - 3) / K_max)) * W_D
    )

    return df


def compute_score(
    df,
    beta=3,
    gamma=1,
    delta=0.8,
):
    N = df["Frequency"].max()
    A_max = df["Abundance"].max()
    A_min = df["Abundance"].min()

    def alpha(A, F):
        return 1 / (1 + np.exp(gamma * (A - F)))

    def score(F, A):
        if F == N:
            return 1
        else:
            a = alpha((A - A_min) / (A_max - A_min), F / N)
            return delta + (1 - delta) * (
                a * (1 - np.exp(-beta * (F / N)))
                + (1 - a) * (1 - np.exp(-beta * (A - A_min) / (A_max - A_min)))
            )

    df["Score"] = df.apply(
        lambda row: score(row["Frequency"], row["Abundance"]), axis=1
    )
    return df


def exponential_decay_penalty(
    df, lambda_cv=10, lambda_rtd=10, cv_thresh=0.3, rtd_thresh=0.05
):
    """
    Calculates the exponential decay penalty for CV and RTD values in a DataFrame.

    This function applies an exponential decay penalty based on predefined thresholds for
    coefficient of variation (CV) and retention time deviation (RTD). Values exceeding
    these thresholds are penalized, reducing their impact in subsequent analysis.

    Args:
        df (pd.DataFrame): DataFrame containing 'CV' and 'RTD' columns.
        lambda_cv (float): Decay rate for CV penalties.
        lambda_rtd (float): Decay rate for RTD penalties.
        cv_thresh (float): Threshold above which CV values are penalized.
        rtd_thresh (float): Threshold above which RTD values are penalized.

    Returns:
        pd.DataFrame: A copy of the input DataFrame with an additional 'Penalty' column,
                      representing the calculated penalty for each row.

    Example:
        >>> data = {'CV': [0.1, 0.2, 0.3], 'RTD': [5, 10, 15]}
        >>> df = pd.DataFrame(data)
        >>> df = exponential_decay_penalty(df, 10, 10, 0.15, 0.05)
    """
    # Calculate penalties for CV and RTD based on the exponential decay formula
    penalty_cv = np.exp(-lambda_cv * np.maximum(0, df["CV"] - cv_thresh))
    penalty_rtd = np.exp(-lambda_rtd * np.maximum(0, df["RTD"] - rtd_thresh))

    # Combine penalties and add to the DataFrame
    df["Penalty"] = penalty_cv * penalty_rtd

    return df


def peptide_confidence_score(Extractedfeatures, CVThreshold, RTDThreshold):
    """
    Calculate the confidence score for each peptide based on its features.

    This function calculates the confidence score for each peptide based on its abundance, frequency, skewness, kurtosis,
    coefficient of variation (CV), and retention time deviation (RTD). It also applies an exponential decay penalty to
    the CV and RTD values based on predefined thresholds.

    Args:
        featureExtractedPeptides (pd.DataFrame): DataFrame containing peptide features.
        CVThreshold (float): Threshold above which CV values are penalized.
        RTDThreshold (float): Threshold above which RTD values are penalized.

    Returns:
        pd.DataFrame: A copy of the input DataFrame with an additional 'ConfidenceScore' column,
                      representing the calculated confidence score for each peptide.

    Example:
        >>> data = {'Abundance': [100, 200, 300], 'Frequency': [0.1, 0.2, 0.3], 'CV': [0.1, 0.2, 0.3], 'RTD': [5, 10, 15]}
        >>> df = pd.DataFrame(data)
        >>> df = peptide_confidence_score(df, 0.15, 0.05)
    """
    # Calculate the confidence score based on peptide features
    Extractedfeatures = calculate_distribution_adjustment_factor(Extractedfeatures)
    Extractedfeatures = compute_score(Extractedfeatures)
    Extractedfeatures = exponential_decay_penalty(
        Extractedfeatures, cv_thresh=CVThreshold, rtd_thresh=RTDThreshold
    )
    Extractedfeatures["ConfidenceScore"] = (
        Extractedfeatures["Score"]
        * Extractedfeatures["D_factor"]
        * Extractedfeatures["Penalty"]
    )
    return Extractedfeatures
