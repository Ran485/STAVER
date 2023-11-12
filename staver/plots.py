import numpy as np
import pandas as pd
import os
import time
import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.express as px
import matplotlib as mpl
import matplotlib.pyplot as plt
import colorsys
import seaborn as sns
import scanpy as sc
import anndata

from random import sample

# python matplotlib export editable PDF
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.dpi"] = 150


def sample_feature(data, sample_rate=0.2, log_transform=False):
    """sample the experiment number from the dataframe"""
    sample_list = sample(set(data.columns), round(len(data.columns) * sample_rate))
    if log_transform:
        data_sample = np.log2(data[sample_list])
    else:
        data_sample = data[sample_list]
    data_sample = data_sample.dropna(axis=0, thresh=1)
    return data_sample


# sample_feature(plasma_data, sample_rate = 0.5, log_transform = True)


def density_plot(
    data=None, curve_type="normal", title=None, x_title=None, y_title=None, outpath=None
):
    """density plot"""

    hist_arr = []
    exp_name = []
    data.columns = data.columns.str.split("_").str[2:5].str.join("_")
    hist_cols = list(data.columns)
    for id, col in enumerate(hist_cols):
        exp_name.append(col)
        col = data.iloc[:, id].values
        col = [x for x in col if not pd.isnull(x)]
        col = np.array(col)
        hist_arr.append(col)
    hist_arr = [i for i in hist_arr]

    group_labels = exp_name
    colors = ["#333F44", "#37AA9C", "#94F3E4"]
    # Create distplot with curve_type set to 'normal'
    fig = ff.create_distplot(
        hist_arr,
        group_labels,
        show_hist=False,
        colors=colors,
        curve_type=curve_type,  #'kde'. 'normal
        show_rug=False,
        histnorm="probability density",
    )
    # Add title information
    # fig.update_layout(title_text='Curve and Rug Plot')
    fig.update_layout(
        title=title,
        xaxis_title=x_title,
        yaxis_title=y_title,
        template="ggplot2",
        font=dict(size=12, color="Black"),
        xaxis=dict(showgrid=True),
        yaxis=dict(showgrid=True),
        plot_bgcolor="#E5ECF6",
    )  # E5ECF6, #fafafa, ededed
    fig.show()
    if outpath is not None:
        fig.write_image(outpath + f"{title}.pdf", engine="kaleido")


def get_sub_string(s):
    parts = s.split("_")
    # 返回第三个（index=2）和第五个（index=4）之间的字符，用'_'连接
    return "_".join(parts[2:5])


def box_plot(
    data,
    log_transform=True,
    id_vars=None,
    x_axis=None,
    y_axis=None,
    title=None,
    x_title=None,
    y_title=None,
    outpath=None,
):
    df_melt = pd.melt(data.reset_index(), id_vars=[id_vars], value_vars=data.columns)
    if log_transform:
        df_melt["value"] = np.log2(df_melt["value"] + 1)
    df_melt[x_axis] = df_melt[x_axis].apply(get_sub_string)
    fig = px.box(df_melt, x=x_axis, y="value")
    # fig.update_layout(title_text='Curve and Rug Plot')
    fig.update_layout(
        title=title,
        xaxis_title=x_title,
        yaxis_title=y_title,
        template="ggplot2",
        font=dict(size=12, color="Black"),
        xaxis=dict(showgrid=True),
        yaxis=dict(showgrid=True),
        plot_bgcolor="#E5ECF6",
    )  # E5ECF6, #fafafa, ededed
    fig.show()
    if outpath is not None:
        fig.write_image(outpath + f"{title}.pdf", engine="kaleido")


def boxplot(
    df, log_transform=True, title=None, x_title=None, y_title=None, outpath=None
):
    fig = go.Figure()
    df_ = df.copy()
    if log_transform:
        df_ = np.log2(df_ + 1)
    df_.columns = df_.columns.str.split("_").str[2:5].str.join("_")
    for col in df_.columns:
        fig.add_trace(
            go.Box(
                y=df_[col],
                name=col,
                line=dict(color="#636EFA"),  # 设置线条颜色
                marker_color="#636EFA",
            )  # 设置点的颜色
        )

    fig.update_layout(
        title=title,
        xaxis_title=x_title,
        yaxis_title=y_title,
        template="ggplot2",
        font=dict(size=12, color="Black"),
        xaxis=dict(showgrid=True),
        yaxis=dict(showgrid=True),
        plot_bgcolor="#E5ECF6",
    )  # E5ECF6, #fafafa, ededed
    fig.show()

    if outpath is not None:
        fig.write_image(outpath + f"{title}.pdf", engine="kaleido")


def generate_color_map(labels, alpha_nodes=1.0, alpha_links=0.55):
    hues = [i / len(labels) for i in range(len(labels))]
    colors_hsv = [(h, 0.5, 0.8) for h in hues]

    colors_rgba_nodes = [colorsys.hsv_to_rgb(*hsv) for hsv in colors_hsv]
    colors_rgba_links = [(r, g, b, alpha_links) for r, g, b in colors_rgba_nodes]
    colors_rgba_nodes = [(r, g, b, alpha_nodes) for r, g, b in colors_rgba_nodes]

    colors_str_nodes = [
        "rgba({},{},{},{})".format(int(r * 255), int(g * 255), int(b * 255), a)
        for r, g, b, a in colors_rgba_nodes
    ]
    colors_str_links = [
        "rgba({},{},{},{})".format(int(r * 255), int(g * 255), int(b * 255), a)
        for r, g, b, a in colors_rgba_links
    ]

    return dict(zip(labels, colors_str_nodes)), dict(zip(labels, colors_str_links))


def draw_sankey(data, filename=None):
    all_labels = sorted(
        list(
            set(
                data["Tissue_type"].unique().tolist()
                + data["Cancer_type"].unique().tolist()
                + data["Cancer_subtype"].unique().tolist()
            )
        )
    )

    node_color_map, link_color_map = generate_color_map(all_labels)

    # Aggregate data
    links = (
        data.groupby(["Tissue_type", "Cancer_type", "Cancer_subtype"])
        .size()
        .reset_index(name="count")
    )

    final_source = []
    final_target = []
    final_value = []

    for index, row in links.iterrows():
        final_source.append(all_labels.index(row["Tissue_type"]))
        final_target.append(all_labels.index(row["Cancer_type"]))
        final_value.append(row["count"])

        final_source.append(all_labels.index(row["Cancer_type"]))
        final_target.append(all_labels.index(row["Cancer_subtype"]))
        final_value.append(row["count"])

    link_colors = [link_color_map[all_labels[s]] for s in final_source]

    fig = go.Figure(
        go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=all_labels,
                color=[node_color_map[label] for label in all_labels],
            ),
            link=dict(
                source=final_source,
                target=final_target,
                value=final_value,
                color=link_colors,
            ),
        )
    )

    fig.update_layout(
        title_text="Sankey Diagram of Cancer Types and Subtypes", font_size=10
    )

    if filename:
        fig.write_image(filename)

    fig.show()


# 使用示例
# df = pd.DataFrame(data)
# draw_sankey(result_df, filename="sankey_diagram-3.pdf")


def correlation_heatmap(
    data, method="pearson", log_transformed=False, outpath=None, filename=None
):
    """Correaltion heatmap and correlation matrix
    Args:
    --------------
    data -> Dataframe: dataframe of raw data
    method -> str: pearson or spearman
    log_transformed -> bool: whether to use log10 transformed data

    Return:
    -----------
    Dataframe: Correaltion of experiments data

    """
    if log_transformed:
        data = np.log10(data)
    corr = data.corr(method=method)
    ## Platelet selection refernce: https://learnku.com/articles/39890
    plt.figure(figsize=(10, 10))
    # optional cmap: RdYlGn; RdYlGn_r; YlGn; rocket_r; YlGnBu; YlGnBu_r; YlOrBr; YlOrBr_r; YlOrRd; YlOrRd_r; RdBu_r
    sns.clustermap(
        data,
        cmap="Reds",
        annot=False,
        robust=False,
        col_cluster=False,
        row_cluster=False,
        linewidths=0.005,
        fmt=".2f",
    )
    # sn.heatmap(corr, annot=True, cmap='vlag')  ## optional cmap: RdYlGn; RdYlGn_r; YlGn; rocket_r; YlGnBu; YlGnBu_r; YlOrBr; YlOrBr_r; YlOrRd; YlOrRd_r;
    if outpath:
        plt.savefig(outpath + f"/{filename}_corr_heatmap.pdf")
        corr.to_csv(outpath + f"/{filename}_corr_matrix.csv")
    plt.show()
    # return corr


def plot_protein_numbers(data, x_col, y_col, title, y_lim=None, height=6, aspect=1.3):
    """
    Creates a combined box and strip plot for protein number visualization.

    Args:
        data (pd.DataFrame): DataFrame containing the data to be plotted.
        x_col (str): Column name in `data` to be plotted on the x-axis.
        y_col (str): Column name in `data` to be plotted on the y-axis.
        title (str): Title of the plot.
        y_lim (tuple, optional): Tuple specifying the limits for the y-axis (e.g., (0, 7000)).
        height (float, optional): Height (in inches) of each facet. Defaults to 6.
        aspect (float, optional): Aspect ratio of each facet, so that aspect * height gives the width of each facet in inches. Defaults to 1.3.

    Example:
        >>> plot_protein_numbers(protein_num, 'Type', 'Identification of Protein numbers',
                                'Protein numbers of HEK293T QC samples', y_lim=(0, 7000), height=5, aspect=1)
    """
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)

    # Create a box plot
    box_plot = sns.catplot(
        x=x_col,
        y=y_col,
        hue=x_col,
        kind="box",
        legend=False,
        height=height,
        aspect=aspect,
        data=data,
    )

    # Overlay with a strip plot
    sns.stripplot(
        x=x_col,
        y=y_col,
        hue=x_col,
        jitter=True,
        dodge=True,
        marker="o",
        palette="Set2",
        alpha=0.5,
        data=data,
    )

    # Set additional plot attributes
    box_plot.fig.suptitle(title)  # Set the title for the figure
    plt.legend(loc="lower right")
    if y_lim:
        plt.ylim(y_lim)

    plt.show()


def plot_molecular_variance(df, column_name):
    """
    Plots the density curves of the original and STAVER processed data for comparison.

    Args:
        df: A pandas dataframe containing the data.
        column_name: A string representing the column name of the data to be plotted.

    Returns:
        None
    """
    plt.figure(figsize=(5, 4))

    # Density plot of the original data
    sns.kdeplot(
        df[df[column_name] == "Raw data"]["Coefficient of Variation [%]"],
        label="Original Density",
        color="blue",
        linestyle="--",
    )

    # Density plot of the STAVER-processed data
    sns.kdeplot(
        df[df[column_name] == "STAVER"]["Coefficient of Variation [%]"],
        label="STAVER Density",
        color="red",
    )

    plt.legend()
    plt.title("Coefficient of Variation density curve")
    plt.xlabel("Coefficient of Variation [%]")
    plt.ylabel("Density")
    plt.show()


def visualize_cancer_subtypes(
    proteomics_data_path,
    subtypes_file_path,
    color_params=["Batch", "Tissue_type", "Cancer_type"],
    figsize_per_plot=(4.5, 4),
    legend_loc="on data",
    outpath="./",
    filename=None,
    add_outline=False,
    show_plot=False,
):
    """
    Visualize the given proteomics data and cancer subtypes using t-SNE and UMAP.

    Args:
        proteomics_data_path (str): Path to the proteomics data (samples as rows, proteins as columns).
        subtypes_file_path (str): Path to the file with samples and their corresponding cancer subtypes.
        color_params (list): List of parameters to color by in the plots.
        figsize_per_plot (tuple): Size of each individual plot.
        outpath (str): Directory to save the plots.
        filename (str): Prefix of the saved files.

    Returns:
        None. Generates visualization images for t-SNE and UMAP.
    """

    # Load data
    data = pd.read_csv(proteomics_data_path, index_col=0).T.replace(np.nan, 0)

    # Load subtype information
    subtypes = pd.read_csv(subtypes_file_path, index_col=0)

    # Check for missing subtype information
    if not all(sample in subtypes.index for sample in data.index):
        raise ValueError("Some samples lack subtype information!")

    # Convert data to AnnData format
    adata = anndata.AnnData(X=data)

    # Add subtype information to AnnData
    adata.obs = subtypes.loc[adata.obs_names]

    # Data normalization (if required)
    sc.pp.scale(adata)

    # Compute t-SNE and UMAP
    sc.tl.tsne(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Set filenames if not provided
    if not filename:
        filename = os.path.basename(proteomics_data_path).split(".")[0]

    # # Visualize t-SNE
    # fig, axes = plt.subplots(figsize=(len(color_params)*figsize_per_plot[0], figsize_per_plot[1]), nrows=1, ncols=len(color_params))
    # for i, param in enumerate(color_params):
    #     sc.pl.tsne(adata, color=param, legend_loc=legend_loc, title=f"t-SNE colored by {param}", ax=axes[i], show=show_plot)
    # fig.tight_layout()
    # fig.savefig(os.path.join(outpath, f"{filename}_tsne.pdf"))
    # plt.close(fig)

    # # Visualize UMAP
    # fig, axes = plt.subplots(figsize=(len(color_params)*figsize_per_plot[0], figsize_per_plot[1]), nrows=1, ncols=len(color_params))
    # for i, param in enumerate(color_params):
    #     sc.pl.umap(adata, color=param, legend_loc=legend_loc, title=f"UMAP colored by {param}", ax=axes[i], show=show_plot)
    # fig.tight_layout()
    # fig.savefig(os.path.join(outpath, f"{filename}_umap.pdf"))
    # plt.close(fig)

    # Determine the number of plots (number of color_params x 2 for both t-SNE and UMAP)
    nrows = 2
    ncols = len(color_params)

    # Create a figure with the necessary number of subplots
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(ncols * figsize_per_plot[0], nrows * figsize_per_plot[1]),
    )

    # Visualize t-SNE on the first row
    for i, param in enumerate(color_params):
        sc.pl.tsne(
            adata,
            color=param,
            legend_loc=legend_loc,
            title=f"t-SNE colored by {param}",
            add_outline=add_outline,
            ax=axes[0, i],
            show=show_plot,
        )

    # Visualize UMAP on the second row
    for i, param in enumerate(color_params):
        sc.pl.umap(
            adata,
            color=param,
            legend_loc=legend_loc,
            title=f"UMAP colored by {param}",
            add_outline=add_outline,
            ax=axes[1, i],
            show=show_plot,
        )

    # Save the combined visualization
    fig.tight_layout()
    fig.savefig(os.path.join(outpath, f"{filename}_combined.pdf"))
    plt.show()
    plt.close(fig)

    return adata


def cancer_specifity_barplot(
    df,
    target_col,
    log_transform=False,
    x_label=None,
    y_label=None,
    title=None,
    save_path=None,
):
    """
    Plots a barplot.

    Args:
        df (pd.DataFrame): The dataframe.
        target_col (str): The column name to be plotted as the target column.
        x_label (str): The x-axis label. Default is None.
        y_label (str): The y-axis label. Default is None.
        title (str): The title of the plot. Default is None.
        save_path (str): The path to save the image. Default is None, indicating no saving.

    Returns:
        None
    """

    plt.figure(figsize=(4, 2.7))
    if log_transform:
        sns.barplot(x=df.index, y=np.log2(df[target_col] + 1))
    # Use seaborn's barplot function to plot the graph
    else:
        sns.barplot(x=df.index, y=df[target_col])

    # Set the title and axis labels
    if title:
        plt.title(title)
    if x_label:
        plt.xlabel(x_label)
    else:
        plt.xlabel("Index")
    if y_label:
        plt.ylabel(y_label)
    else:
        plt.ylabel(target_col)

    plt.xticks(rotation=30, ha="right")  # Rotate the x-axis labels for better display

    # Save the image (if save_path is specified)
    if save_path:
        plt.savefig(save_path, bbox_inches="tight")

    # Show the image
    plt.show()


def plot_combined_barplots(
    datasets, target_col, titles, log_transforms, x_labels, y_labels, save_path=None
):
    """
    Plots multiple barplots on a single figure, one for each dataset provided.

    This function creates a series of bar plots for different datasets, allowing for
    comparisons across them. It supports log transformation and allows customization
    of plot titles, x and y labels.

    Args:
        datasets (list of pd.DataFrame): List of DataFrames to plot.
        target_col (str): Name of the column in DataFrame to plot.
        titles (list of str): Titles for each subplot.
        log_transforms (list of bool): Whether to apply log transformation for each dataset.
        x_labels (list of str): X-axis labels for each subplot.
        y_labels (list of str): Y-axis labels for each subplot.
        save_path (str, optional): Path to save the figure. If None, the figure is not saved.

    Raises:
        ValueError: If the length of any list argument does not match the length of 'datasets'.

    Examples:
        >>> dataset1 = pd.DataFrame({'CTSE': [1, 2, 3]})
        >>> dataset2 = pd.DataFrame({'CTSE': [4, 5, 6]})
        >>> plot_combined_barplots(
                datasets=[dataset1, dataset2],
                target_col="CTSE",
                titles=["Dataset 1", "Dataset 2"],
                log_transforms=[True, False],
                x_labels=["Index", "Index"],
                y_labels=["Value", "Value"]
            )
    """

    if not all(
        len(lst) == len(datasets)
        for lst in [titles, log_transforms, x_labels, y_labels]
    ):
        raise ValueError(
            "All list arguments must be of the same length as the 'datasets' list."
        )

    n = len(datasets)
    fig, axes = plt.subplots(1, n, figsize=(n * 4, 4.2))

    if n == 1:  # Handle the case of a single dataset
        axes = [axes]

    for i, (df, title, log_transform, x_label, y_label) in enumerate(
        zip(datasets, titles, log_transforms, x_labels, y_labels)
    ):
        ax = axes[i]
        plot_data = np.log2(df[target_col] + 1) if log_transform else df[target_col]
        sns.barplot(ax=ax, x=df.index, y=plot_data)
        ax.set_title(title)
        ax.set_xlabel(x_label if x_label else "Index")
        ax.set_ylabel(y_label if y_label else target_col)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")

    plt.subplots_adjust(left=0.1)  # Adjust the left margin

    if save_path:
        plt.savefig(save_path, bbox_inches="tight")

    plt.tight_layout()
    plt.show()


def compare_execution_time(functions, data, labels=None):
    """
    Measures the execution time of multiple functions and plots a bar chart.

    Args:
        functions (list): A list of functions to be tested. These functions should accept the same arguments.
        data: The input data to be passed to each function. Ensure the data format is suitable for all functions.
        labels (list of str, optional): A list of labels for each function. Defaults to the names of the functions.

    Raises:
        ValueError: If the number of provided labels does not match the number of functions.

    Usage:
        >>> def test_func1(data):
        ...     # Simulate some computation
        ...     time.sleep(1)
        >>> def test_func2(data):
        ...     # Simulate some other computation
        ...     time.sleep(2)
        >>> data = [1, 2, 3]  # Data for testing the functions
        >>> compare_execution_time([test_func1, test_func2], data, ['Function 1', 'Function 2'])

    This function will output a bar chart of the execution times for test_func1 and test_func2.
    """
    if labels is None:
        labels = [f.__name__ for f in functions]
    elif len(labels) != len(functions):
        raise ValueError("Number of labels must match number of functions.")

    times = []
    for func in functions:
        start_time = time.time()
        if func.__name__ == "maxLFQ":
            func(data)
        else:
            func(data.to_numpy())
        elapsed_time = time.time() - start_time
        times.append(elapsed_time)

    plt.bar(labels, times, color=["blue", "green", "red", "purple", "orange"])
    plt.ylabel("Execution Time (seconds)")
    plt.title("Execution Time Comparison")
    plt.show()


def plot_performance(measure_func, sample_sizes):
    """
    Measures and plots the performance of a function over varying sample sizes.

    This function measures the execution time of a specified function for different sample sizes and
    plots the performance curve.

    Args:
        measure_func (callable): A function that measures the performance. It should accept a list of sample sizes as input.
        sample_sizes (list of int): A list of sample sizes to measure the performance over.

    The function plots the execution time against the sample sizes using matplotlib.

    Usage:
        >>> def measure_performance(sizes):
        ...     # Function to measure performance for different sample sizes
        ...     # Return a list of execution times
        ...     pass
        >>> sample_sizes = [10, 50, 100, 500, 1000]
        >>> plot_performance(measure_performance, sample_sizes)
    """
    times = measure_func(sample_sizes)

    plt.plot(sample_sizes, times, marker="o")
    plt.xlabel("Sample Size")
    plt.ylabel("Time (seconds)")
    plt.title("Function Performance with Varying Sample Sizes")
    plt.grid(True)
    plt.show()
