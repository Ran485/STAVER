#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : protein_dynamic_range.py
@Time     : 2023/05/02 00:08:22
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 21112030023@m.fudan.edu.cn
@License  : (C)Copyright 2022-2023, DingLab-CHINA-SHNAGHAI
@Function : None
'''
# here put the import lib


import csv
import os
import time

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42


def intensity_plot(data, project, color, highlight_proteins, axis_max):
    """
    This function plots the intensity of proteins
    Parameters:

        data (dataframe): dataframe containing protein intensities
        project (string): name of the project
        color (string): color of the plot
        highlight_proteins (list): list of proteins to be highlighted
        axis_max (int): maximum value of x-axis

    Returns:
        None
    """
    protein_dict = {}
    mean_intensity = np.log10(data.mean(axis=1).sort_values(ascending=False)).to_frame()
    mean_intensity.rename(
        columns={0: 'Quantative protein intensity[Log10]'}, inplace=True
    )
    mean_intensity['Rank of quantified proteins'] = (
        np.arange(0, len(mean_intensity.index)) + 1
    )
    mean_intensity['color'] = color
    for protein in highlight_proteins:
        try:
            a = mean_intensity.loc[protein]
            mean_intensity.loc[protein, 'color'] = 'blue'
            protein_list = mean_intensity.loc[
                protein,
                ['Rank of quantified proteins', 'Quantative protein intensity[Log10]'],
            ].to_list()
            protein_dict[protein] = tuple(protein_list)
        except Exception:
            print(protein + ' is not exsited in file!')

    plt.rcParams["figure.dpi"] = 400

    ax = mean_intensity.plot.scatter(
        figsize=(7, 5),
        x='Rank of quantified proteins',
        y='Quantative protein intensity[Log10]',
        s=1,
        c=mean_intensity['color'],
    )

    protein_plot = mean_intensity[mean_intensity['color'] == 'blue']

    protein_plot.plot.scatter(
        x='Rank of quantified proteins',
        y='Quantative protein intensity[Log10]',
        s=2,
        c=protein_plot['color'],
        ax=ax,
    )

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.title(project, fontsize=10)
    plt.xlim(xmax=axis_max)
    plt.xlabel('Rank of quantified proteins', fontsize=10)
    plt.ylabel('Quantative protein intensity[Log10]', fontsize=10)

    for k, v in protein_dict.items():
        new_v = tuple(x + 0.15 for x in list(v))
        ax.annotate(str(k), xy=v, xytext=new_v, color='k', size=5)
    # plt.show()
    formats = ['.png', '.pdf', '.eps']
    # formats = ['.pdf']
    for format_N in formats:
        plt.savefig(
            output + 'proteins_' + project + '_' + time_name + format_N,
            transparent=True,
        )


# proteins=['ALB','APOA1','HPX','APOE','C1QC','VWF','GOT1','HLA-DRA','IL1R2','SEMA3F','TNNC1','MDH2','GLG1']
# Representative plasma proteins in plasma across different orders of magnitude
represent_proteins = [
    'ALB',
    'APOA1',
    'FGB',
    'AFM',
    'C2',
    'F12',
    'ADIPOQ',
    'CRP',
    'SEPP1',
    'LBP',
    'CEACAM1',
    'HPX',
    'APOE',
    'C1QC',
    'VWF',
    'GOT1',
    'PARK7',
    'ACE',
    'IL1RAP',
    'HLA-DRA',
    'IL1R2',
    'SEMA3F',
    'TNNI3',
    'PNLIP',
    'AFP',
    'VEGFA',
    'IL1B',
    'IL6',
    'EPCAM',
    'PSA',
    'AFP',
    'EGFR',
    'VEGF',
    'CSF1',
    'CEA',
    'CSF1',
]

# Drug targets that have been reported in the literature
drudged_proteins = [
    'ERBB2',
    'ERBB3',
    'ERBB4',
    'EGFR',
    'KITLG',
    'PDGFRB',
    'PDGFRA',
    'FLT4',
    'FGFR3',
    'FGFR2',
    'FGFR4',
    'FGFR1',
    'ACVR1B',
    'BMPR1A',
    'ACVR1',
    'MAP3K1',
    'ROS1',
    'TNFRSF1A',
    'MTOR',
    'MRC1',
    'ITGAV',
    'IGF1R',
    'ABCC1',
    'MMP15',
    'MMP14',
    'TLR2',
    'CD44',
    'ITGA2B',
    'IL1R2',
    'PDGFRB',
    'MMP16',
    'SLC12A2',
    'EPCAM',
    'ICAM1',
]


if __name__ == '__main__':
    # 使用有空值有表达量的csv文件
    # inpath = r'/Volumes/T7_Shield/staver/results/DIA_repeat20/processed_data/Protein_Dynamic.csv'
    inpath = r'/Users/ranpeng/Library/Mobile Documents/com~apple~CloudDocs/Desktop/CTC项目/processed_data/protein_matrix_dynamic.csv'
    # The results of the analysis are saved in the following path
    output = r'/Users/ranpeng/Library/Mobile Documents/com~apple~CloudDocs/Desktop/CTC项目/figs/'
    if not os.path.isdir(output):
        os.mkdir(output)
    time_name = time.strftime('%Y%m%d', time.localtime(time.time()))

    data = pd.read_csv(inpath, index_col=[0], header=[0, 1])
    data.columns.names = [None, None]
    data.index = data.index.get_level_values(0)

    # data.drop(['Gene_Id','P_value','mean0','mean1','Fold_change'],axis=1,inplace=True)

    Cohort1 = 'Exp147561'
    Cohort2 = 'Exp147560'
    # Cohort3 = 'Cohort3'

    intensity_plot(data, 'All', '#82D2F4', drudged_proteins, axis_max=800)
    intensity_plot(data[Cohort1], Cohort1, '#B9EF84', drudged_proteins, axis_max=800)
    intensity_plot(data[Cohort2], Cohort2, '#F94E07', drudged_proteins, axis_max=800)
