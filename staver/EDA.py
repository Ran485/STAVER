import pandas as pd
import numpy as np
import missingno as msno
import seaborn as sns
import matplotlib.pyplot as plt
from IPython.display import display

"""
Module for exploratory data analysis (EDA).

This module provides functions for checking missing values, data types, and visualizing numerical and categorical features.

Functions:
    check_missing_value(data): Check the distribution of missing values in a dataset.
    check_dtypes_missing(data): Check the data types and missing values in a DataFrame.
    get_missing_rate(data): Get the missing rate of features in a dataset.
    get_missing_le50(data): Get the features with missing rates less than or equal to 50%.
    split_data_dtypes(data): Split the features into numerical and categorical columns.
    split_numerical_serial_fea(data, features): Split the numerical features into continuous and discrete variables.
    Visualization_numerical_distribution(data, numerical_serial_fea, col_wrap=5): Visualize the distribution of numerical continuous features.
    plot_category_features(data, category_columns, r, c): Plot the distribution of categorical features.

File: /Users/ranpeng/Desktop/github/STAVER/staver/EDA.py
"""


# ------ Check missing values and distributions in the dataset -----------
def check_missing_value(data):
    """
    Check the distribution of missing values in a dataset.

    Args:
        data: The DataFrame containing the data.
    """
    # Check the total number of rows with missing data
    print(
        f"There are {data.isnull().any().sum()} columns in the dataset with missing values."
    )
    # Visualization of missing values
    msno.matrix(data)


def check_dtypes_missing(data):
    """
    Check the data types and missing values in a DataFrame.

    This function provides information about the column types and the number of missing values in the DataFrame.

    Args:
        data: The DataFrame to be checked.

    Returns:
        None
    """
    tab_info = pd.DataFrame(data.dtypes).T.rename(index={0: "Data Type"})
    tab_info = tab_info.append(
        pd.DataFrame(data.isnull().sum()).T.rename(index={0: "Missing Values"})
    )
    tab_info = tab_info.append(
        pd.DataFrame(data.isnull().sum() / data.shape[0] * 100).T.rename(
            index={0: "Missing Rate (%)"}
        )
    )
    print(
        "-" * 10
        + " Displaying information about column types and missing values "
        + "-" * 10
    )
    display(tab_info)


def get_missing_rate(data):
    """
    Get the missing rate of features in a dataset.

    Args:
        data: The DataFrame containing the data.

    Returns:
        DataFrame: The missing rate of features.
    """
    missing = data.isnull().sum() / len(data)
    missing = missing[missing > 0]
    missing.sort_values(inplace=True)
    missing.plot.bar()

    # Visualization of missing rate
    df = pd.DataFrame(missing).rename(columns={"index": "col", 0: "missing_rate"})
    # sns.kdeplot(df['missing_rate'],shade=True) # density plot
    return df


def get_missing_le50(data):
    """
    Get the features with missing rates less than or equal to 50%.

    Args:
        data: The DataFrame containing the data.

    Returns:
        DataFrame: The features with missing rates less than or equal to 50%.
    """
    feature_have_null_dict = (data.isnull().sum() / len(data)).to_dict()
    feature_null_lessThanHalf = []
    for key, value in feature_have_null_dict.items():
        if value < 0.5:
            feature_null_lessThanHalf.append(key)
    null_lessThanHalf_data = data.loc[:, feature_null_lessThanHalf]
    print(
        "The columns of features with less than 50 percent missing values: {}".format(
            null_lessThanHalf_data.columns
        )
    )
    return null_lessThanHalf_data


# ------ Check data types: 1) Numerical features 2) Categorical features -----------
# - Features are generally composed of categorical and numerical features, and numerical
# features can be further divided into continuous and discrete variables.
# - Categorical features sometimes have non-numeric relationships and sometimes have numeric
# relationships. For example, in 'grade', the grades A, B, C, etc., whether they are just
# categories or if A is better than others, need to be analyzed and judged.
# - Numerical features can be directly used for modeling, but often need to be binned,
# transformed into WOE encoding, and then used for standard scorecard operations. From the
# perspective of model performance, feature binning is mainly to reduce the complexity of
# variables, reduce the impact of variable noise on the model, and improve the correlation
# between independent variables and dependent variables, making the model more stable.


def split_data_dtypes(data):
    """
    Split the features into numerical and categorical columns.

    Args:
        data: The DataFrame containing the data.

    Returns:
        tuple: A tuple containing the lists of numerical and categorical columns.
    """
    numerical_columns = list(data.select_dtypes(exclude=["object"]).columns)
    category_columns = list(
        filter(lambda x: x not in numerical_columns, list(data.columns))
    )
    print(f"There are {len(numerical_columns)} numerical columns in the dataset.")
    print(f"There are {len(category_columns)} category columns in the dataset.")
    print("The list of numerical columns are: {}".format(numerical_columns))
    print("The list of category columns are: {}".format(category_columns))
    return numerical_columns, category_columns


def split_numerical_serial_fea(data, features):
    """
    Split the numerical features into continuous and discrete variables.

    Args:
        data: The DataFrame containing the data.
        features: The list of numerical features.

    Returns:
        tuple: A tuple containing the lists of continuous and discrete numerical features.
    """
    numerical_serial_fea = []
    numerical_noserial_fea = []
    for fea in features:
        temp = data[fea].nunique()
        if temp <= 10:
            numerical_noserial_fea.append(fea)
            continue
        numerical_serial_fea.append(fea)
    print(
        f"There are {numerical_serial_fea} numerical serial features in the {data} dataset."
    )
    print(
        f"There are {numerical_noserial_fea} numerical non-serial features in the {data} dataset."
    )
    return numerical_serial_fea, numerical_noserial_fea


# - Check the distribution of a specific numerical feature, and see if the variable follows a normal distribution.
# If the variable does not follow a normal distribution, it can be log-transformed to see if it conforms to a normal distribution.
# - If you want to standardize a batch of data, you must exclude the data that has already been normalized.
# - The reason for normalization: In some cases, normality or non-normality can make the model converge faster.
# Some models require data to be normal (e.g., GMM, KNN). It is important to ensure that the data is not too skewed,
# as extreme skewness may affect the model's predictive results.
def Visualization_numerical_distribution(data, numerical_serial_fea, col_wrap=5):
    """
    Visualize the distribution of numerical continuous features.

    Args:
        data: The DataFrame containing the data.
        numerical_serial_fea: The list of continuous numerical features.
        col_wrap: The number of columns to wrap the plots.

    Returns:
        None
    """
    f = pd.melt(data, value_vars=numerical_serial_fea)
    g = sns.FacetGrid(f, col="variable", col_wrap=col_wrap, sharex=False, sharey=False)
    g = g.map(sns.distplot, "value")


def plot_category_features(data, category_columns, r, c):
    """
    Plot the distribution of categorical features.

    Args:
        data: The DataFrame containing the data.
        category_columns: The list of categorical columns.
        r: The number of rows in the plot grid.
        c: The number of columns in the plot grid.

    Returns:
        None
    """
    plt.figure()
    plt.figure(figsize=(15, 20))
    i = 1
    for col in category_columns:
        # r represents the number of rows, c represents the number of columns, so there are a total of r*c plots
        plot_envs = plt.subplot(r, c, i)
        i += 1
        v = data[col].value_counts().reset_index()[:10]
        fig = sns.barplot(x=v["index"], y=v[col])
        for item in fig.get_xticklabels():
            item.set_rotation(45)
        plt.title(col)
    plt.tight_layout()
    # plt.show()
