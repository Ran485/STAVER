import pandas as pd
import numpy as np


def max_lfq_quantities(
    peptide_data: pd.DataFrame,
    protein_column: str = 'Protein',
    peptide_col: str = 'Peptide',
) -> pd.DataFrame:
    """
    Calculate MaxLFQ quantities for the given peptide data.

    Parameters:
    -----------
    peptide_data : pandas.DataFrame
        A DataFrame containing peptides as rows and samples as columns.
        Each cell should contain the intensity of the peptide in the corresponding sample.
    protein_column : str, optional
        The name of the column in peptide_data that contains the corresponding protein for each peptide. (default is 'Protein')

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing MaxLFQ quantities for each protein in each sample.
    """
    # Step 1: Log2 transformation of the intensities
    sample_columns = peptide_data.columns[peptide_data.columns != protein_column]
    sample_columns = sample_columns[sample_columns != peptide_col]
    log2_peptide_data = np.log2(peptide_data[sample_columns])

    # Step 2: Calculate the median intensity of each peptide across all samples
    median_peptide_intensity = log2_peptide_data.median(axis=1)

    # Step 3: Calculate the relative intensity of each peptide in each sample
    relative_peptide_intensity = log2_peptide_data.sub(median_peptide_intensity, axis=0)

    # Step 4: Group the peptides by protein and calculate the median relative intensity for each protein in each sample
    max_lfq_quantities = (
        relative_peptide_intensity.join(peptide_data[protein_column])
        .groupby(protein_column)
        .median()
    )

    # Step 5: Convert the MaxLFQ quantities back to linear space
    max_lfq_quantities = 2**max_lfq_quantities

    return max_lfq_quantities

# In this optimized version, we use Pandas' built-in functions, such as join and groupby, 
# n order to take full advantage of Pandas' performance. Also.
# We add an optional parameter called protein_column so that the user can specify the name
# of the column containing the protein information. This optimized version should
# be able to compute MaxLFQ volumes faster when working with large datasets.
if __name__ == '__main__':
    # Load the peptide data
    protein_data = pd.read_csv("/Volumes/T7_Shield/staver/test_data/test.csv")

    # Calculate the MaxLFQ quantities
    max_lfq_quantities = max_lfq_quantities(protein_data)

    # Save the MaxLFQ quantities to a CSV file
    # max_lfq_quantities.to_csv('data/max_lfq_quantities.csv')
