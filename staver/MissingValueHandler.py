#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : MissingValueHandler.py
@Time     : 2023/05/02 00:07:44
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
'''
# here put the import lib


import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer, KNNImputer
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler


class MissingValueHandler:
    """
    A class to handle missing values in proteomics data.

    Attributes:
    -----------
    data : pandas.DataFrame
        The input data containing missing values.

    Methods:
    --------
    drop_rows_with_missing_values(threshold: float = None) -> pandas.DataFrame
        Removes rows containing any missing values or based on a threshold.

    impute_missing_values(strategy: str = 'mean') -> pandas.DataFrame
        Imputes missing values using SimpleImputer with specified strategy.

    knn_impute_missing_values(n_neighbors: int = 5) -> pandas.DataFrame
        Imputes missing values using KNNImputer with specified number of neighbors.

    max_value_imputation() -> pandas.DataFrame
        Imputes missing values using the maximum value in each column.

    min_value_imputation() -> pandas.DataFrame
        Imputes missing values using the minimum value in each column.

    distribution_imputation(std_multiplier: float = 1) -> pandas.DataFrame
        Imputes missing values based on the distribution of the data.
    """

    @staticmethod
    def remove_rows_below_threshold(
        data: pd.DataFrame, threshold: float = 0.3
    ) -> pd.DataFrame:
        """
        Removes rows containing missing values above the given threshold.

        Parameters:
        -----------
        threshold : float, optional
            The percentage threshold of non-missing values required in a row to be kept. (default is 0.8)

        Returns:
        --------
        pandas.DataFrame
            Data with rows containing missing values above the threshold removed.
        """
        return data.dropna(thresh=threshold * len(data.columns))

    @staticmethod
    def simple_impute_values(
        data: pd.DataFrame, strategy: str = "mean"
    ) -> pd.DataFrame:
        """
        Imputes missing values using SimpleImputer with specified strategy.

        Parameters:
        -----------
        strategy : str, optional
            The imputation strategy, can be 'mean', 'median', or 'most_frequent'. (default is 'mean')

        Returns:
        --------
        pandas.DataFrame
            Imputed data with missing values filled.
        """
        imputer = SimpleImputer(missing_values=np.nan, strategy=strategy)
        imputed_data = imputer.fit_transform(data)
        return pd.DataFrame(imputed_data, columns=data.columns)

    @staticmethod
    def knn_impute_missing_values(
        data: pd.DataFrame, n_neighbors: int = 5
    ) -> pd.DataFrame:
        """
        Imputes missing values using KNNImputer with specified number of neighbors.

        Parameters:
        -----------
        n_neighbors : int, optional
            Number of neighboring samples to use for imputing missing values. (default is 5)

        Returns:
        --------
        pandas.DataFrame
            Imputed data with missing values filled.
        """
        imputer = KNNImputer(n_neighbors=n_neighbors)
        imputed_data = imputer.fit_transform(data)
        return pd.DataFrame(imputed_data, columns=data.columns, index=data.index)

    @staticmethod
    def max_impute_values(data: pd.DataFrame) -> pd.DataFrame:
        """Imputes missing values using the maximum value in each column."""
        max_values = data.max()
        imputed_data = data.fillna(max_values)
        return imputed_data

    @staticmethod
    def min_impute_values(data: pd.DataFrame) -> pd.DataFrame:
        """Imputes missing values using the minimum value in each column."""
        # One tenth of the minimum value
        min_values = data.min() / 10
        imputed_data = data.fillna(min_values)
        return imputed_data

    @staticmethod
    def distribution_impute_values(
        data: pd.DataFrame, std_multiplier: float = 1
    ) -> pd.DataFrame:
        """
        Imputes missing values based on the distribution of the data.

        Parameters:
        -----------
        std_multiplier : float, optional
            Standard deviation multiplier for generating imputed values. (default is 1)

        Returns:
        --------
        pandas.DataFrame
            Imputed data with missing values filled.
        """

        def fillna_with_distribution(
            col: pd.Series, std_multiplier: float
        ) -> pd.Series:
            mean = col.mean()
            std = col.std()
            null_count = col.isnull().sum()
            if null_count > 0:
                col.fillna(
                    value=np.random.normal(
                        loc=mean, scale=std * std_multiplier, size=null_count
                    ),
                    inplace=True,
                )
            return col

        imputed_data = data.apply(fillna_with_distribution, args=(std_multiplier,))
        return imputed_data

    @staticmethod
    def impute_missing_values(data: pd.DataFrame) -> pd.DataFrame:
        """
        Imputes missing values in a pandas DataFrame using scikit-learn's SimpleImputer and KNNImputer classes.

        Args:
            data (pd.DataFrame): The DataFrame to impute missing values in.

        Returns:
            pd.DataFrame: The imputed DataFrame.
        """
        # Define the imputation methods
        numerical_imputer = SimpleImputer(strategy="mean")
        categorical_imputer = KNNImputer(n_neighbors=5)

        # Define the column transformer
        transformer = ColumnTransformer(
            transformers=[
                (
                    "num",
                    numerical_imputer,
                    data.select_dtypes(include="number").columns,
                ),
                (
                    "cat",
                    categorical_imputer,
                    data.select_dtypes(include="object").columns,
                ),
            ]
        )

        # Define the pipeline
        pipeline = Pipeline([("imputer", transformer), ("scaler", StandardScaler())])
        # Fit the pipeline to the data
        pipeline.fit(data)
        # Transform the data
        imputed_data = pipeline.transform(data)
        # Convert the transformed data back to a DataFrame
        imputed_data = pd.DataFrame(imputed_data, columns=data.columns)

        return imputed_data
