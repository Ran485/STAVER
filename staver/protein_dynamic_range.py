#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : protein_dynamic_range.py
@Time     : 2022/11/10 19:08:54
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
import json
import time
import os
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42

def intensity_plot(data, project, color, highlight_proteins, axis_max, output, time_name):
    """
    Plots the intensity of proteins as a scatter plot.

    Args:
        data (pd.DataFrame): DataFrame containing protein intensities.
        project (str): Name of the project for the plot title.
        color (str): Color for the plot points.
        highlight_proteins (list of str): List of proteins to be highlighted in the plot.
        axis_max (int): Maximum value for the x-axis.
        output (str): Output directory to save the plot.

    Returns:
        None
    """
    mean_intensity = np.log10(data.mean(axis=1).sort_values(ascending=False)).to_frame()
    mean_intensity.rename(columns={0: 'Quantitative protein intensity[Log10]'}, inplace=True)
    mean_intensity['Rank of quantified proteins'] = np.arange(len(mean_intensity)) + 1
    mean_intensity['color'] = color

    for protein in highlight_proteins:
        if protein in mean_intensity.index:
            mean_intensity.loc[protein, 'color'] = 'blue'

    ax = mean_intensity.plot.scatter(
        x='Rank of quantified proteins', y='Quantitative protein intensity[Log10]',
        s=1, c=mean_intensity['color'], figsize=(7, 5)
    )

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.title(project, fontsize=10)
    plt.xlim(xmax=axis_max)
    plt.xlabel('Rank of quantified proteins', fontsize=10)
    plt.ylabel('Quantitative protein intensity[Log10]', fontsize=10)

    for protein in highlight_proteins:
        if protein in mean_intensity.index:
            position = mean_intensity.loc[protein]
            ax.annotate(protein, xy=(position['Rank of quantified proteins'], position['Quantitative protein intensity[Log10]']),
                        xytext=(5, 0), textcoords='offset points', ha='right', va='bottom', fontsize=5)

    for fmt in ['.png', '.pdf', '.eps']:
        plt.savefig(os.path.join(output, f'proteins_{project}_{time_name}{fmt}'), transparent=True)
    plt.show()


def plot_protein_intensity(inpath, outpath, config_path=None):
    """
    Process and plot protein intensity data.

    Args:
        inpath (str): Input file path for protein data.
        outpath (str): Output directory path for plots.
        config_path (str): Path to the protein configuration file.

    Returns:
        None
    """
    os.makedirs(outpath, exist_ok=True)
    time_name = time.strftime('%Y%m%d', time.localtime())

    protein_lists = protein_list_config(config_path)
    data = pd.read_csv(inpath, index_col=0, header=[0, 1])
    data.replace(np.nan, data.min(), inplace=True)
    data.columns.names = [None, None]
    data.index = data.index.get_level_values(0)

    for cohort, color in [('All', '#82D2F4'), ('Good_Circulation', '#B9EF84'), ('Bad_Circulation', '#F94E07')]:
        intensity_plot(data.get(cohort, data), cohort, color, protein_lists['represent_proteins'], 6500, outpath, time_name)


def protein_list_config(inpath = None):
    """
    Loads protein lists from a JSON configuration file.

    Args:
        config_path (str): Path to the JSON configuration file.

    Returns:
        dict: Dictionary containing protein lists.
    """
    # protein_config.json
    if inpath:
        with open(inpath, 'r') as f:
            return json.load(f)
    else:
        return {
            # Representative plasma proteins in plasma across different orders of magnitude 
            "represent_proteins": [
                "ALB", "APOA1", "FGB", "AFM", "C2", "F12", "ADIPOQ", "CRP", "SEPP1", "LBP",
                "CEACAM1", "HPX", "APOE", "C1QC", "VWF", "GOT1", "PARK7", "ACE", "IL1RAP",
                "HLA-DRA", "IL1R2", "SEMA3F", "TNNI3", "PNLIP", "AFP", "VEGFA", "IL1B", "IL6",
                "EPCAM", "PSA", "AFP", "EGFR", "VEGF", "CSF1", "CEA", "CSF1"
            ],
            # Drug targets that have been reported in the literature
            "drugged_proteins": [
                "ERBB2", "ERBB3", "ERBB4", "EGFR", "KITLG", "PDGFRB", "PDGFRA", "FLT4", "FGFR3",
                "FGFR2", "FGFR4", "FGFR1", "ACVR1B", "BMPR1A", "ACVR1", "MAP3K1", "ROS1", "TNFRSF1A",
                "MTOR", "MRC1", "ITGAV", "IGF1R", "ABCC1", "MMP15", "MMP14", "TLR2", "CD44", "ITGA2B",
                "IL1R2", "PDGFRB", "MMP16", "SLC12A2", "EPCAM", "ICAM1"
            ]
        }

# Example usage
if __name__ == '__main__':
    inpath = '/path/to/Protein_matrix_rawdata_sorted_dynamic.csv'
    outpath = '/path/to/output/figs/'
    
    # inpath = r'/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/validation_data/merge-data/Protein_matrix_rawdata_sorted_dynamic.csv'
    # # The results of the analysis are saved in the following path
    # outpath = r'/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/validation_data/figs-1/'
    plot_protein_intensity(inpath, outpath)
