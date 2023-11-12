#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : protein_accumulation_plots.py
@Time     : 2022/11/10 19:08:34
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
'''
# here put the import lib

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
import os
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)


def accumulate_data(accu_data):
    """
    Calculate cumulative count of unique non-zero values in each column.

    Args:
        data (pd.DataFrame): DataFrame with data to be accumulated.

    Returns:
        pd.DataFrame: DataFrame with accumulated data.
    """
    arr = np.array(accu_data.T.values)
    a = arr.cumsum(0)
    data = [sum(a[0] == 1)]
    for i in range(1, (len(a))):
        cz = arr[int(i)] - a[int(i) - 1]
        data.append(sum(cz == 1))
    data = np.array(data)
    data_accu = data.cumsum()
    df_accu = pd.DataFrame(data_accu, index=accu_data.columns)
    df_accu = df_accu.reset_index()
    df_accu.rename(columns={0: 'accumulation'}, inplace=True)
    return df_accu


def count_non_null(data):
    """
    Count non-null values in each column of a DataFrame.

    Args:
        data (pd.DataFrame): DataFrame with data to be counted.

    Returns:
        pd.DataFrame: DataFrame with counts of non-null values.
    """
    data_count = data.count()
    data_count = data_count.reset_index()
    data_count.rename(columns={0: 'count'}, inplace=True)
    return data_count


def plot_number_of_proteins(count, accu, colors, outpath, time_name):
    """
    Plot number of quantified proteins and cumulative protein numbers.

    Args:
        count (pd.DataFrame): DataFrame containing count data.
        accu (pd.DataFrame): DataFrame containing accumulation data.
        colors (dict): Dictionary of colors for each cohort.
        outpath (str): Output path for saving the plot.
        time_name (str): Time identifier for the output file.
    """
    fig, ax1 = plt.subplots(figsize=(14, 9))
    sns.barplot(x='Exp', y='count', data=count, hue='cohort', palette=colors, ax=ax1)
    ax1.set_xlabel('Cohort', fontsize=14)
    ax1.set_ylabel('Quantified protein numbers', fontsize=14)
    ax1.legend(loc='upper left', fontsize=10, shadow=True, fancybox=True)
    plt.ylim(0, count['count'].max() * 1.1)

    ax2 = ax1.twinx()
    # sns.lineplot(x='Exp', y='accumulation', data=accu, palette=colors, marker='o', linewidth=2, ax=ax2, label='Your Label')
    sns.lineplot(x='Exp', y='accumulation', data=accu, palette=colors, marker='o', linewidth=2, ax=ax2)
    ax2.set_ylabel('Cumulative protein numbers', fontsize=14)
    ax2.legend().remove()

    plt.title('Quantified protein numbers', fontsize=14)
    plt.savefig(f'{outpath}/GPs_numbers_{time_name}.pdf', dpi=400, transparent=True)
    plt.show()


def plot_boxenplot(data, colors, outpath, time_name):
    """
    Plot a boxenplot for the given data.

    Args:
        data (pd.DataFrame): DataFrame containing the data to plot.
        colors (dict): Dictionary of colors for each cohort.
        outpath (str): Output path for saving the plot.
        time_name (str): Time identifier for the output file.
    """
    fig, ax = plt.subplots(figsize=(21, 8))
    sns.boxenplot(x='Exp', y='iFOT', hue='cohort', data=data, palette=colors, ax=ax)
    ax.set_xlabel('Sample', fontsize=14)
    ax.set_ylabel('LFQ intensity[log10]', fontsize=14)
    ax.legend().remove()
    plt.savefig(f'{outpath}/boxplot_{time_name}.pdf')
    plt.show()


def process_and_plot_data(inpath, outpath, time_name, colors):
    """
    Processes proteomics data and generates plots.

    Args:
        inpath (str): Path to the input CSV file.
        outpath (str): Path to the output directory for saving plots.
        time_name (str): A time-based string to append to the output filenames.
        colors (dict): Dictionary of colors for the plot.

    Returns:
        None
    """
    os.makedirs(outpath, exist_ok=True)

    data = pd.read_csv(inpath, header=[0, 1], index_col=[0, 1])
    accu_data = (data > 0).astype(int)

    data_count = count_non_null(data)
    data_accumulation = accumulate_data(accu_data)

    plot_number_of_proteins(data_count, data_accumulation, colors, outpath, time_name)

    # Prepare data for boxenplot
    gene_ID = data.index.get_level_values(1)
    boxplot_data = pd.DataFrame(data.values, index=gene_ID, columns=data.columns)
    boxplot_data = boxplot_data.unstack().reset_index()
    boxplot_data.rename(columns={0: 'iFOT'}, inplace=True)
    boxplot_data['iFOT'] = np.log10(boxplot_data['iFOT'])

    plot_boxenplot(boxplot_data, colors, outpath, time_name)


if __name__ == '__main__':
    time_name = time.strftime('%Y%m%d', time.localtime(time.time()))
    # inpath = '/path/to/matrix_data_FOT_sorted_dynamic.csv'
    # outpath = '/path/to/output/sfig-1-1'
    inpath = r'/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/validation_data/merge-data/Protein_matrix_rawdata_sorted_dynamic.csv'
    # The results of the analysis are saved in the following path
    outpath = r'/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/validation_data/figs-1/'
    colors = {
        'Good_Circulation': '#D53B60',
        'Bad_Circulation': '#77236D',
    }

    process_and_plot_data(inpath, outpath, time_name, colors)
