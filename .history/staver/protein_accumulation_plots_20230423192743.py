import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
import os

# python matplotlib PDF 不断字
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42


# 计算累计值
def Accumulate(accu_data):
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


def Count(data):
    data_count = data.count()
    data_count = data_count.reset_index()
    data_count.rename(columns={0: 'count'}, inplace=True)
    return data_count


def Number_plot(count, accu):
    fig = plt.figure(figsize=(14, 9))
    ax1 = sns.barplot(x='Exp', y='count', hue='cohort', data=count, palette=colors)
    ax1.set_xlabel('Cohort', fontsize=14)
    ax1.set_ylabel('Quantified protein numbers', fontsize=14)
    ax1.legend(loc='upper left', fontsize=10, shadow=True, fancybox=True)
    ax1_max = count['count'].max() * 1.1
    plt.ylim(0, ax1_max)

    # ax2 = ax1.twinx()
    # ax2 = sns.lineplot(x='Exp', y='accumulation', hue='cohort', data=accu, marker='o', linewidth=2, palette=colors)
    # ax2.set_ylabel('Cumulative protein numbers', fontsize=14)
    # ax2.legend().remove()
    # ax2_max = accu['accumulation'].max() * 1.1
    # plt.ylim(0, ax2_max)
    accu.to_csv(outpath + '/accumulation_protein_numbers.csv')

    ax2 = ax1.twinx()
    ax2 = sns.lineplot(
        x='Exp',
        y='accumulation',
        hue='cohort',
        data=accu,
        palette=colors,
        marker='o',
        drawstyle='default',
        linewidth=2,
    )
    ax2.set_ylabel('Cumulative protein numbers', fontsize=14)
    ax2.legend().remove()
    plt.xticks([])

    plt.xticks([])
    plt.title('Quantified protein numbers', fontsize=14)  # Protein identification
    # for ax in [ax1, ax2]:
    #     for orient in ['top', 'left', 'right']:
    #         ax.spines[orient].set_visible(False)
    plt.savefig(
        outpath + '/GPs_numbers_' + time_name + '.pdf', dpi=400, transparent=True
    )
    plt.show()


def Boxnplot(plot_data):
    fig = plt.figure(figsize=(21, 9))
    ax = sns.boxenplot(x='Exp', y='iFOT', hue='cohort', data=plot_data, palette=colors)
    plt.xticks([])
    ax.legend().remove()
    ax.set_xlabel('Sample', fontsize=14)
    ax.set_ylabel('LFQ intensity[log10]', fontsize=14)
    plt.savefig(outpath + 'boxplot_' + time_name + '.pdf')
    plt.show()


colors = {
    # 'Normal': '#77236D',
    # 'Covid_19': '#D53B60',
    # 'H7N9': '#F97C6F'
    'STAVER': '#D53B60',  # 77236D
    'Raw Data': '#77236D',
    # 'H7N9': '#F97C6F'
}

if __name__ == '__main__':
    time_name = time.strftime('%Y%m%d', time.localtime(time.time()))
    inpath = r'/Volumes/T7_Shield/staver/results/DIA_repeat20/processed_data/Protein_Dynamic.csv'
    outpath = r'/Volumes/T7_Shield/staver/results/DIA_repeat20/processed_data/figs-2//'
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    data = pd.read_csv(inpath, header=[0, 1], index_col=[0, 1])
    data[data <= 0] = np.nan

    accu_data = data.copy()
    accu_data.replace(np.nan, 0, inplace=True)
    accu_data[accu_data > 0] = 1

    data_count = Count(data)
    data_accumulation = Accumulate(accu_data)
    Number_plot(data_count, data_accumulation)

    # boxplot
    gene_ID = data.index.get_level_values(1)
    boxplot_data = pd.DataFrame(data.values, index=gene_ID, columns=data.columns)
    boxplot_data = boxplot_data.unstack()
    bp = pd.DataFrame(boxplot_data.values, index=boxplot_data.index)
    bp = bp.reset_index()
    bp.rename(columns={0: 'iFOT'}, inplace=True)
    bp['iFOT'] = np.log10(bp['iFOT'])
    Boxnplot(bp)
