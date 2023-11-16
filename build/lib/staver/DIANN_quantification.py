import os
from collections import defaultdict

import pandas as pd
import numpy as np
import time


def IbaqData(ibaq_file):
    global ibaq_data, gene_data
    ibaq_data = pd.read_csv(ibaq_file)
    gene_data = ibaq_data.drop("Proteinid", axis=1)[~(ibaq_data["Symbol"].isnull())]
    gene_data.drop_duplicates("Symbol", keep="first", inplace=True)
    gene_data.drop("Cleave_Count", axis=1, inplace=True)


def FindFiles(tsv_path):
    csvs = []
    for root, dirs, files in os.walk(tsv_path, topdown=False):
        csvs.extend(
            os.path.join(root, names)
            for names in files
            if names.endswith(".tsv")
            and "gene" not in names
            and "stats" not in names
            and names.startswith(".") == False
        )
    return csvs


def TsvMerge(csvs):
    tsv_merge = pd.DataFrame(
        [],
        columns=["Protein.Ids", "Precursor.Normalised", "Stripped.Sequence", "Q.Value"],
    )
    for gene_file in csvs:
        data = pd.read_csv(
            gene_file,
            sep="\t",
            usecols=[
                "Protein.Ids",
                "Precursor.Normalised",
                "Q.Value",
                "Stripped.Sequence",
                "Precursor.Id",
            ],
        )
        data = data[data["Q.Value"] <= q]
        # data.drop('Q.Value',axis=1,inplace=True)
        data = data[~(data["Protein.Ids"].isnull())]
        data["LibraryCount"] = 1
        tsv_merge = pd.concat([tsv_merge, data], axis=0, sort=False)
    tsv_merge = tsv_merge.groupby(
        ["Protein.Ids", "Stripped.Sequence", "Precursor.Id"]
    ).agg({"Precursor.Normalised": np.mean, "LibraryCount": np.sum})
    tsv_merge = tsv_merge[tsv_merge["LibraryCount"] >= 2]  # 选择被大于1个库鉴定到的PSM

    tsv_merge = pd.DataFrame(tsv_merge)
    tsv_merge.reset_index(inplace=True)
    tsv_merge.drop("Precursor.Id", axis=1, inplace=True)

    return tsv_merge


def TsvHandle(csvs):
    data = TsvMerge(csvs)
    data = data[data["Protein.Ids"].map(lambda x: not x.startswith("rev_"))]
    data["PSMs"] = 1
    # data = data[data['Q.Value']<= q]
    protein_ids_count = (
        data["Protein.Ids"]
        .str.split(";", expand=True)
        .count(axis=1)
        .rename("Protein.Ids.Count")
    )
    data_split_row = (
        data.drop("Protein.Ids", axis=1)
        .join(
            data["Protein.Ids"]
            .str.split(";", expand=True)
            .stack()
            .reset_index(level=1, drop=True)
            .rename("Protein.Ids")
        )
        .join(protein_ids_count)
    )
    data_split_row["Area"] = (
        data_split_row["Precursor.Normalised"] / data_split_row["Protein.Ids.Count"]
    )
    data_split_row["Unique Peptide Num"] = data_split_row["Protein.Ids.Count"]
    data_split_row.loc[
        data_split_row["Unique Peptide Num"] != 1, "Unique Peptide Num"
    ] = 0
    # data_split_row['Unique Peptide Num'][data_split_row['Unique Peptide Num']!=1]=0

    peptide_grouped = data_split_row.groupby(["Stripped.Sequence", "Protein.Ids"])[
        ["Area", "PSMs", "Unique Peptide Num"]
    ].sum()
    peptide_grouped.loc[
        peptide_grouped["Unique Peptide Num"] != 0, "Unique Peptide Num"
    ] = 1
    # peptide_grouped['Unique Peptide Num'][peptide_grouped['Unique Peptide Num']!=0]=1
    peptide_grouped["Peptide Num"] = 1
    protein_grouped = peptide_grouped.groupby("Protein.Ids")[
        ["Area", "PSMs", "Unique Peptide Num", "Peptide Num"]
    ].sum()

    protein_grouped["Protein Num"] = 1

    data_merge = pd.merge(
        protein_grouped,
        ibaq_data[["Proteinid", "Symbol", "Cleave_Count"]],
        how="inner",
        left_index=True,
        right_on="Proteinid",
    )

    data_merge["Proteinid"] = data_merge["Proteinid"].apply(lambda x: f"{str(x)};")
    # print(data_merge)
    # data_merge.rename(columns=lambda x:x.upper(), inplace=True)
    # print(data_merge.columns)
    # columns_list=list(set(data_merge.columns)-{'Cleave_count','Protein.Ids','PG.Normalised','gene_protein_count','Sequence'})
    agg_dict = {k: sum for k in list(set(data_merge.columns) - {"Symbol"})}

    gene_grouped = data_merge.groupby("Symbol").agg(agg_dict)

    gene_grouped["iBAQ"] = gene_grouped["Area"] / gene_grouped["Cleave_Count"]
    gene_grouped["FoT(1e-5)"] = gene_grouped["iBAQ"] * 1e5 / gene_grouped["iBAQ"].sum()
    gene_grouped.drop("Cleave_Count", axis=1, inplace=True)
    return pd.merge(
        gene_grouped,
        gene_data,
        how="inner",
        left_index=True,
        right_on="Symbol",
    )


def tsv(path, out_path):
    folders = [
        path + f
        for f in os.listdir(path)
        if os.path.isdir(path + f) and f != "quantification"
    ]
    for tsv_path in folders:
        csvs = FindFiles(tsv_path)
        name = tsv_path.split("/")[-1].split(".")[0]
        if len(csvs) != 327:
            print(name + " has not enough files!")
            continue
        elif os.path.isfile(out_path + name + ".csv"):
            print(name + " has been merged!")
            continue
        else:
            print(name)
            file_out = TsvHandle(csvs)
            file_out.to_csv(out_path + name + ".csv", index=False, columns=columns)


if __name__ == "__main__":
    import sys

    pd.set_option("display.max_columns", None)
    time_name = time.strftime("%Y%m%d", time.localtime(time.time()))

    columns = [
        "Entrezgene",
        "Symbol",
        "Name",
        "Proteinid",
        "Area",
        "iBAQ",
        "FoT(1e-5)",
        "Protein Num",
        "Peptide Num",
        "Unique Peptide Num",
        "PSMs",
        "Hgnc",
        "Map_Location",
        "Pfam",
        "Taxid",
        "Type_Of_Gene",
        "Unigene",
        "Genomic_Pos_Hg19.Chr",
        "Genomic_Pos_Hg19.Start",
        "Genomic_Pos_Hg19.End",
        "Genomic_Pos_Hg19.Strand",
        "Pharos.Target_Id",
        "Umls.Cui",
        "Ec",
    ]
    q = 0.01
    path = (
        "/Volumes/Disk3/STAVER-revised/validation20_293T/validation20_293T_327_FDR001/"
    )
    out_path = "/Volumes/Disk3/STAVER-revised/validation20_293T/merge-data/"
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    # path = sys.argv[1]
    # out_path = sys.argv[2]
    if not os.path.isdir(out_path):
        os.makedirs(out_path)
    print(f"The input path is {path}")
    print(f"The output path is {out_path}")
    ibaq_file = r"/Volumes/T7_Shield/DIA-QC/results/paper_clinical_infomation//Uniprot_9606_iBAQ_20200203.csv"
    IbaqData(ibaq_file)
    tsv(path, out_path)
