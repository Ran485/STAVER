import pandas as pd
import numpy as np
import networkx as nx


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


def peptide_ratio(peptide_intensity1, peptide_intensity2):
    temp_df = pd.DataFrame({'intensity1': peptide_intensity1, 'intensity2': peptide_intensity2})
    return temp_df.apply(lambda row: row['intensity1'] / row['intensity2'] if row['intensity1'] > 0 and row['intensity2'] > 0 else (np.nan if row['intensity1'] > 0 or row['intensity2'] > 0 else 0), axis=1)

def protein_ratio_matrix(protein_data):
    sample_names = protein_data.columns[2:]
    sample_number = len(sample_names)
    ratio_matrix = np.zeros((sample_number, sample_number))

    for i in range(sample_number):
        for j in range(i, sample_number):
            peptide_intensity1 = protein_data[sample_names[i]]
            peptide_intensity2 = protein_data[sample_names[j]]
            peptide_ratios = peptide_ratio(peptide_intensity1, peptide_intensity2)
            max_peptide_ratio = np.nanmax(peptide_ratios)
            if not np.isnan(max_peptide_ratio):
                ratio_matrix[i, j] = max_peptide_ratio
                ratio_matrix[j, i] = 1 / max_peptide_ratio
    return ratio_matrix

def compute_relative_abundance(tree, root, sample_number):
    abundance = np.ones(sample_number)
    visited = set()
    stack = [(root, 1)]

    while stack:
        node, weight = stack.pop()
        visited.add(node)
        abundance[node] = weight

        for neighbor in tree.neighbors(node):
            if neighbor not in visited:
                edge_weight = tree[node][neighbor]['weight']
                stack.append((neighbor, weight * edge_weight))

    return abundance

def protein_abundance(protein_data):
    protein_names = protein_data["Protein"].unique()
    sample_names = protein_data.columns[2:]
    sample_number = len(sample_names)
    abundance_matrix = np.zeros((len(protein_names), sample_number + 2))

    for i, protein_name in enumerate(protein_names):
        protein_subdata = protein_data[protein_data["Protein"] == protein_name]
        ratio_matrix = protein_ratio_matrix(protein_subdata)
        graph = nx.from_numpy_matrix(ratio_matrix)
        tree = nx.minimum_spanning_tree(graph)
        root = np.argmax(np.sum(ratio_matrix, axis=1))
        abundance = compute_relative_abundance(tree, root, sample_number)
        total_intensity = protein_subdata[sample_names].sum().sum()
        normalized_abundance = abundance * total_intensity / abundance.sum()
        abundance_matrix[i, 2:] = normalized_abundance

    abundance_df = pd.DataFrame(abundance_matrix, index=protein_names, columns=['Protein', 'Peptide'] + sample_names.tolist())
    abundance_df['Protein'] = protein_names
    abundance_df['Peptide'] = protein_data.groupby('Protein')['Peptide'].apply(lambda x: ','.join(x))

    return abundance_df

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


    # # Load the test data and compute protein abundances
    # test_data = pd.read_csv("test_data.csv")
    # result = protein_abundance(test_data)