import pandas as pd
import numpy as np


def normalize_data(data):
    """
    Normalize data to range between 0 and 1.

    Args:
        data (pd.Series or np.array): Data to be normalized.

    Returns:
        pd.Series or np.array: Normalized data.
    """
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def validate_data(data):
    """
    Validates data for numeric values and handles missing values.

    Args:
        data (pd.DataFrame): Data to be validated.

    Raises:
        ValueError: If data is empty or contains non-numeric values.

    Returns:
        pd.DataFrame: Validated and potentially modified DataFrame.
    """
    if data is None or data.empty:
        raise ValueError("Data is empty, please check the input file.")

    if not data.applymap(np.isreal).all().all():
        raise ValueError(
            "Data contains non-numeric values, please check the input file."
        )

    if data.isnull().any().any():
        print("Warning: Null value detected, already filled with zeros.")
        data = data.fillna(0)

    return data


def sigmoid(x):
    """
    Computes the sigmoid function.

    Args:
        x (float or np.array): Input value(s).

    Returns:
        float or np.array: Sigmoid function output.
    """
    return 1 / (1 + np.exp(-x))


def peptide_scoring(
    data, weight_variation=1, weight_frequency=5, CV_thresh=0.5, non_linear=False
):
    """
    Calculate peptide scores based on frequency and variation.

    Args:
        data (pd.DataFrame): Input data containing peptide information.
        weight_variation (int): Weight factor for variation. Defaults to 1.
        weight_frequency (int): Weight factor for frequency. Defaults to 5.
        CV_thresh (int): Coefficient of variation threshold. Defaults to 0.5.
        non_linear (bool): Apply sigmoid function for non-linear transformation. Defaults to False.

    Returns:
        pd.DataFrame: DataFrame with calculated peptide scores.
    """
    data = validate_data(data)
    peptide_frequency = data.count(axis=1) / len(data.columns)
    peptide_variation = data.std(axis=1) / data.mean(axis=1)
    normalized_frequency = normalize_data(peptide_frequency)
    normalized_variation = normalize_data(peptide_variation)

    peptide_scores = (
        weight_frequency
        * normalized_frequency
        * (1 - weight_variation * normalized_variation)
    )
    peptide_scores *= (peptide_variation > 0) & (peptide_frequency > 0)
    peptide_scores += (peptide_variation == 0) & (peptide_frequency > 0)
    penalty = (peptide_variation > CV_thresh) & (data.count(axis=1) <= 3)
    peptide_scores *= 1 - penalty

    if non_linear:
        peptide_scores = sigmoid(peptide_scores)

    peptide_scores = normalize_data(peptide_scores)

    return pd.DataFrame(
        {
            "Peptide": data.index,
            "Count": data.count(axis=1),
            "Frequency": peptide_frequency,
            "Variation": peptide_variation,
            "Score": peptide_scores,
        }
    ).sort_values(by="Score", ascending=False)


def filter_data(data, sample_threshold=0.1, variation_threshold=0.2):
    """
    Filter peptide data based on frequency and variation thresholds.

    Args:
        data (pd.DataFrame): Input data containing peptide information.
        sample_threshold (float): Frequency threshold for filtering. Defaults to 0.1.
        variation_threshold (float): Variation threshold for filtering. Defaults to 0.2.

    Returns:
        pd.DataFrame: Filtered peptide data.
    """
    return data[
        (data["Frequency"] >= sample_threshold)
        & (data["Variation"] <= variation_threshold)
    ]


if __name__ == "__main__":
    data_path = "/path/to/Protein_matrix_1.csv"
    output_path = "/path/to/Protein_matrix_score_result.csv"
    data = pd.read_csv(data_path, index_col=0)
    scored_data = peptide_scoring(data, non_linear=True)
    scored_data.to_csv(output_path, index=False)
