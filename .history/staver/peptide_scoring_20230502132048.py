#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : peptide_scoring.py
@Time     : 2023/05/02 00:07:54
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
'''
# here put the import lib


import pandas as pd
import numpy as np


def normalize_data(data):
    """Normalize data to range between 0 and 1"""
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def validate_data(data):
    """Data validation and error correction mechanisms"""
    if data is None or data.empty:
        raise ValueError("Data is empty, please check the input file.")

    if not data.applymap(np.isreal).all().all():
        raise ValueError(
            "Data contains non-numeric values, please check the input file."
        )

    if data.isnull().any().any():
        # data = data.fillna(0)
        print("Warning: Null value detected, already filled with zeros.")

    return data


def sigmoid(x):
    """
    Sigmoid 函数

    Parameters:
    -----------
    x : float or np.array
        输入值

    Returns:
    --------
    float or np.array
        Sigmoid 函数输出值
    """
    return 1 / (1 + np.exp(-x))


def peptide_scoring(
    data, weight_variation=1, weight_frequency=1, CV_thresh=2, non_linear=False
):
    """
    Calculate peptide scores.

    Parameters:
    -----------
    data : pandas.DataFrame
        The input data containing peptide information.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the peptide scores.

    Raises:
    -------
    ValueError
        If the input data is empty or contains non-numeric values.
    """
    # 验证数据
    data = validate_data(data)
    # 计算肽段在样本中出现的频率
    peptide_count = data.count(axis=1)
    peptide_frequency = data.count(axis=1) / len(data.columns)
    # 计算每个肽段丰度的变异
    peptide_variation = data.std(axis=1) / data.mean(axis=1)
    # 对鉴定频次和变异进行归一化处理
    normalized_frequency = normalize_data(peptide_frequency)
    normalized_variation = normalize_data(peptide_variation)

    # 以变异值为0，鉴定频次为所有样本的肽段的最高得分设定为1，其他的则小于1
    # peptide_scores = normalized_frequency * (1 - normalized_variation)
    # 考虑变异权重更大，肽段鉴定频次影响较小
    peptide_scores = (
        weight_frequency
        * normalized_frequency
        * (1 - weight_variation * normalized_variation)
    )

    # Set variation to 0 and identification frequency to 1 for all peptides
    peptide_scores *= (peptide_variation > 0) & (peptide_frequency > 0)
    peptide_scores += (peptide_variation == 0) & (peptide_frequency > 0)

    # Penalize peptides with high variation and low identification frequency
    penalty = (peptide_variation > CV_thresh) & (peptide_count <= 3)
    peptide_scores *= 1 - penalty
    if non_linear:
        # 使用sigmoid函数进行非线性变换
        peptide_scores = sigmoid(peptide_scores)

    # 限定评分范围在 [0-1] 区间内
    peptide_scores = normalize_data(peptide_scores)

    # 保存结果
    result = pd.DataFrame(
        {
            'Peptide': data.index,
            'Count': peptide_count,
            'Frequency': peptide_frequency,
            'Variation': peptide_variation,
            'Score': peptide_scores,
        }
    )

    return result.sort_values(by='Score', ascending=False)


def filter_data(data, sample_threshold=0.1, variation_threshold=0.2):
    """
    Filter peptide data based on frequency and variation thresholds.

    Parameters:
    -----------
    data : pandas.DataFrame
        The input data containing peptide information.
    sample_threshold : float, optional
        The frequency threshold for filtering peptides (default is 0.5).
    variation_threshold : float, optional
        The variation threshold for filtering peptides (default is 0.1).

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the filtered peptide data.
    """
    # 根据频率阈值筛选肽段
    filtered_peptides = data[data['Frequency'] >= sample_threshold]
    # 根据变异阈值筛选肽段
    filtered_peptides[data['Variation'] <= variation_threshold]

    return filtered_peptides


if __name__ == '__main__':
    # 使用函数计算结果并保存
    data = pd.read_csv(
        '/Volumes/T7_Shield/staver/results/DIA_repeat20/processed_data/Protein_matrix_1.csv',
        index_col=0,
    )
    result = peptide_scoring(data, non_linear=True)
    result.to_csv(
        '/Volumes/T7_Shield/staver/results/DIA_repeat20/processed_data/Protein_matrix_score_result.csv',
        index=False,
    )
