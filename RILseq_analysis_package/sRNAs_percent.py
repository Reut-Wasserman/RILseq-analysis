from RILseq_analysis_package.utils import *
import argparse
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import matplotlib
import os


def find_highest_sRNAs(df, sRNAs_amount):
    """
    Returns a list of the sRNAs with the highest percentages.
    :param df: a data frame of the percentages of each sRNA in each condition.
    """
    max_values = df.max().sort_values(ascending=False)
    return max_values.head(sRNAs_amount).index.tolist()


def get_counts_table(sRNAs, path, experiments, replicates):
    """
    Returns a data frame of the number of reads of each sRNA for each condition. The number is the sum of the read
    numbers at the different replicates.
    :param sRNAs: a list of sRNAs names.
    :param experiment: RILseq or RNAseq
    """
    df = pd.read_csv(path, sep="\t")
    names_dic = {i:i.replace("_cutadapt_bwa", "") for i in df.columns}
    names_dic["Unnamed: 0"] = "gene"
    df = df.rename(columns=names_dic)

    df = df[df["gene"].isin(sRNAs)]
    df.index = df["gene"].values
    relevant_cols = []
    for condition in experiments:
        df[condition] = df[[condition + i for i in replicates]].sum(axis=1)
        relevant_cols.append(condition)

    return df[relevant_cols]


def convert_to_percent(df):
    row_sum = df.sum(axis=1)
    return df.div(row_sum/100, axis=0)


def rename(rna):
    lst = list(rna)
    lst[0] = rna[0].upper()
    return "".join(lst)


def get_colors(genes_list, genes_colors_dic):
    colors = plt.cm.tab20c.colors + plt.cm.tab20b.colors
    res = []
    for gene in genes_list:
        if gene in genes_colors_dic.keys():
            res.append(genes_colors_dic[gene])
        else:
            color = colors[len(genes_colors_dic)]
            res.append(color)
            genes_colors_dic[gene] = color
    return res, genes_colors_dic


def chimeras_bar_plot(sRNAs_amount, sRNAs, experiments, replicates, base_path, only_between_chr_chimeras, rna_seq):
    """
    Creates bar plots of the percentages of the sRNAs for each condition according to the number of chimers, chimeric fragments,
    RNAseq reads and RILSeq reads. The 15 sRNAs with the highest percentages are colored, the other are gray.
    :param sRNAs: a list of the sRNAs names
    :param only_between_chr_chimeras: if True, create the plot according to the between-chromosome chimeras.
    """
    experiments = [i for i in experiments if "wt" not in i]
    df_chimeras_amount = pd.DataFrame(np.zeros((len(experiments), len(sRNAs))), columns=sRNAs, index=experiments)
    df_fragments_amount = pd.DataFrame(np.zeros((len(experiments), len(sRNAs))), columns=sRNAs, index=experiments)
    df_reads_amount_RILseq = pd.DataFrame(np.zeros((len(experiments), len(sRNAs))), columns=sRNAs, index=experiments)
    counts_table_RILseq = get_counts_table(sRNAs, os.path.join(base_path, "all_counts_table.txt"), experiments, replicates)
    if rna_seq is not None:
        df_reads_amount_RNAseq = pd.DataFrame(np.zeros((len(experiments), len(sRNAs))), columns=sRNAs, index=experiments)
        counts_table_RNAseq = get_counts_table(sRNAs, rna_seq["all_counts_table_file"], experiments, replicates)
    for experiment in experiments:
        chimeras_df = pd.read_excel(os.path.join(base_path, r"RILseq_unified_results.xlsx"), sheet_name=experiment)
        if only_between_chr_chimeras:
            chimeras_df = chimeras_df[chimeras_df["RNA1_chromosome"] != chimeras_df["RNA2_chromosome"]]
            # chromosomes = list(chr_dic.keys())
            # chimeras_df = chimeras_df[((chimeras_df["RNA1 chromosome"] == chromosomes[0]) & (chimeras_df["RNA2 chromosome"] == chromosomes[1])) |
            #                           ((chimeras_df["RNA1 chromosome"] == chromosomes[1]) & (chimeras_df["RNA2 chromosome"] == chromosomes[0]))]
        for rna in sRNAs:
            srna_chimeras = chimeras_df[(chimeras_df["RNA1 name"] == rna) | (chimeras_df["RNA2 name"] == rna)]
            df_chimeras_amount.at[experiment, rna] = srna_chimeras.shape[0]
            df_fragments_amount.at[experiment, rna] = srna_chimeras["interactions"].sum()
            df_reads_amount_RILseq.at[experiment, rna] = counts_table_RILseq[experiment][rna]
            if rna_seq is not None:
                df_reads_amount_RNAseq.at[experiment, rna] = counts_table_RNAseq[experiment][rna]
    dic = {"chimeras":df_chimeras_amount, "fragments":df_fragments_amount, "RILseq reads":df_reads_amount_RILseq}
    if rna_seq is not None:
        dic["RNAseq reads"] = df_reads_amount_RNAseq
    genes_colors_dic = {"other sRNAs":"gray"}
    for name, df in dic.items():
        df = convert_to_percent(df)
        fig = plt.figure()
        relevant_sRNAs = find_highest_sRNAs(df, sRNAs_amount)
        other_sRNAs = [i for i in sRNAs if i not in relevant_sRNAs]
        df["other sRNAs"] = df[other_sRNAs].sum(axis=1)
        # df = df.rename(index=experiments_dic)
        colors, genes_colors_dic = get_colors(relevant_sRNAs + ["other sRNAs"], genes_colors_dic)
        # if only_between_chr_chimeras and EXPERIMENT == "lambda":
        #     if name in ("chimeras", "fragments"):
        #         df = df.loc[["infected 30", "infected 60"]]
        ax = df[relevant_sRNAs + ["other sRNAs"]].plot.bar(stacked=True, color=colors)
        handles, labels = ax.get_legend_handles_labels()
        labels = [rename(i) for i in labels]
        plt.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5))
        plt.xticks(rotation=0)
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x))))
        plt.ylabel(f"{name} %")
        fig_name = f"sRNAs_{name}_percent"
        if only_between_chr_chimeras:
            fig_name = f"{fig_name}_only_between_chr_chimeras"
        plt.savefig(os.path.join(base_path, f"{fig_name}.png"), bbox_inches='tight')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("sRNAs_amount")
    parser.add_argument("--RNAseq_counts", default=None)

    args = parser.parse_args()
    config = load_config(args.config)
    sRNAs = get_RNA_types("sRNA", config["rna_types_excel"])
    all_genes = get_annotation(config["annotation_path"], separate_id_name=True)["name"].values.tolist()
    sRNAs = [i for i in sRNAs if i in all_genes]
    chimeras_bar_plot(args.sRNAs_amount, sRNAs, config["experiments"], config["replicates"], config["base_path"], True, args.RNAseq_counts)
    chimeras_bar_plot(args.sRNAs_amount, sRNAs, config["experiments"], config["replicates"], config["base_path"], False, args.RNAseq_counts)





