import pandas as pd
import numpy as np


def max_lfq_quantities(peptide_data: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate MaxLFQ quantities for the given peptide data.

    Parameters:
    -----------
    peptide_data : pandas.DataFrame
        A DataFrame containing peptides as rows and samples as columns.
        Each cell should contain the intensity of the peptide in the corresponding sample.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing MaxLFQ quantities for each protein in each sample.
    """
    # Step 1: Log2 transformation of the intensities
    log2_peptide_data = np.log2(peptide_data)

    # Step 2: Calculate the median intensity of each peptide across all samples
    median_peptide_intensity = log2_peptide_data.median(axis=1)

    # Step 3: Calculate the relative intensity of each peptide in each sample
    relative_peptide_intensity = log2_peptide_data.sub(median_peptide_intensity, axis=0)

    # Step 4: Group the peptides by protein and calculate the median relative intensity for each protein in each sample
    # Note: You need a column in the DataFrame indicating the corresponding protein for each peptide
    max_lfq_quantities = relative_peptide_intensity.groupby(
        peptide_data['Protein']
    ).median()

    # Step 5: Convert the MaxLFQ quantities back to linear space
    max_lfq_quantities = 2**max_lfq_quantities

    return max_lfq_quantities


# 这个max_lfq_quantities函数接受一个包含肽段数据的DataFrame作为输入，其中行表示肽段，列表示样本。
# 该函数首先将肽段数据的强度值进行log2变换，然后计算每个肽段在所有样本中的中位强度。接下来，计算每个
# 肽段在每个样本中的相对强度。最后，根据蛋白质对肽段进行分组，并计算每个蛋白质在每个样本中的最大LFQ量。
# 这些量最后转换回线性空间，以获得每个蛋白质在每个样本中的MaxLFQ量。


import pandas as pd
import numpy as np


def max_lfq_quantities_optimized(
    peptide_data: pd.DataFrame, protein_column: str = 'Protein'
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


def max_lfq_quantities_optimized(
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


# 在这个优化版本中，我们使用Pandas的内置函数，如join和groupby，以便充分利用Pandas的性能优势。同时，
# 我们添加了一个名为protein_column的可选参数，以便用户可以指定包含蛋白质信息的列名。这个优化版本应该
# 能够在处理大型数据集时更快地计算MaxLFQ量。
if __name__ == '__main__':
    # Load the peptide data
    protein_data = pd.read_csv("/Volumes/T7_Shield/staver/test_data/test.csv")

    # Calculate the MaxLFQ quantities
    max_lfq_quantities = max_lfq_quantities_optimized(protein_data)

    # Save the MaxLFQ quantities to a CSV file
    # max_lfq_quantities.to_csv('data/max_lfq_quantities.csv')
