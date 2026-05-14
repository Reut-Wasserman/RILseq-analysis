import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from RILseq_analysis_package.utils import get_annotation, get_RNA_types, load_config
import argparse



def get_gene_annotation(gene, genes, tRNAs, sRNAs, oRNAs):
    if gene in tRNAs:
        return "tRNA"
    if gene in sRNAs:
        return "sRNA"
    if gene in oRNAs:
        return "oRNA"
    if "3UTR" in gene:
        return "3UTR"
    if "5UTR" in gene:
        return "5UTR"
    if "AS" in gene:
        return "AS"
    if "IGR" in gene:
        return "IGR"
    if "IGT" in gene:
        return "IGT"
    if gene in genes:
        return "CDS"
    return "unknown"


def get_RNAs_types(rna_types_excel):
    sRNAs = set(get_RNA_types("sRNA", rna_types_excel))
    tRNAs = set(get_RNA_types("tRNA", rna_types_excel))
    oRNAs = set(get_RNA_types("oRNA", rna_types_excel))
    return sRNAs, tRNAs, oRNAs


def RNA1_RNA2_annotations(chr_dic, base_path, experiments, rna_types_excel, annotation_path, chimeras_or_fragments, chrom=None):
    """
    Creates heatmaps of the chimeras annotations, with and without the annotations numbers.
    :param chimeras_or_fragments: if "chimeras" create the heatmaps according to the amount of chimeras. Else, according
    to the amount of fragments.
    :param chrom: if None, create the heatmaps base on the chimeras between the two chromosomes. Else, create the heatmaps
     base on the chimeras of the given chromosome.
    """
    genes = get_annotation(annotation_path, chromosome=chrom, separate_id_name=True)["name"].values.tolist()
    sRNAs, tRNAs, oRNAs = get_RNAs_types(rna_types_excel)
    if chrom:
        dir_ = chr_dic[chrom] + "annotation"
    else:
        dir_ = "_".join(chr_dic.keys()) + "annotation"
    output_path = os.path.join(base_path, dir_)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    for experiment in experiments:
        # if chrom == CHR2 and treatment == "Hfq_":
        #     continue
        df = pd.read_excel(os.path.join(base_path, "RILseq_unified_results.xlsx"), sheet_name=experiment)
        if chrom is not None:
            df = df[(df["RNA1 chromosome"] == chrom) & (df["RNA2 chromosome"] == chrom)]
        else:
            chr1, chr2 = chr_dic.keys()
            df = df[((df["RNA1 chromosome"] == chr1) & (df["RNA2 chromosome"] == chr2)) |
                    ((df["RNA1 chromosome"] == chr2) & (df["RNA2 chromosome"] == chr1))]
        df["RNA1_annotation"] = df["RNA1 name"].apply(get_gene_annotation, genes=genes, tRNAs=tRNAs, sRNAs=sRNAs, oRNAs=oRNAs)
        df["RNA2_annotation"] = df["RNA2 name"].apply(get_gene_annotation, genes=genes, tRNAs=tRNAs, sRNAs=sRNAs, oRNAs=oRNAs)
        cur_options = list(set(df["RNA1_annotation"].values.tolist() + df["RNA2_annotation"].values.tolist()))
        options = ["3UTR", "5UTR", "CDS", "AS", "IGR", "sRNA", "tRNA", "IGT"] #"oRNA"

        if chimeras_or_fragments == "chimeras":
            df["i"] = df.index
            df = df[["RNA1_annotation", "RNA2_annotation", "i"]]
            small_df = df.groupby(["RNA1_annotation", "RNA2_annotation"]).count()
            small_df_pairs = small_df.index
            amount = small_df["i"].values
        else:
            df = df[["RNA1_annotation", "RNA2_annotation", "interactions"]]
            small_df = df.groupby(["RNA1_annotation", "RNA2_annotation"]).sum()
            small_df_pairs = small_df.index
            amount = small_df["interactions"].values
        dic = {}
        for i in range(len(amount)):
            dic[small_df_pairs[i]] = amount[i]
        results = pd.DataFrame(0, index=options, columns=options)
        for option1 in options:
            for option2 in options:
                if (option1, option2) in small_df_pairs:
                    results[option1][option2] = dic[(option1, option2)]
        for annot in [True, False]:
            plt.figure()
            sns.heatmap(results, annot=annot, cmap='Blues', fmt=".6g", annot_kws={"size": 8})
            plt.xlabel("RNA1")
            plt.ylabel("RNA2")
            with_annot = "with_annot" if annot else "without_annot"
            plt.savefig(os.path.join(output_path, f"{experiment}_{chimeras_or_fragments}_{with_annot}.png"))


def get_RNAs_type_fractions(df, options, chimeras_or_fragments):
    """
    :param df: chimeras data-frame
    :param options: the annotations
    :param chimeras_or_fragments: if "chimeras" calculate the fraction according to the amount of chimeras. Else, according
    to the amount of fragments.
    :return: a list of the fraction of each annotation
    """
    amount_dic = {i:0 for i in options}
    for RNA in "12":
        small_df = df[[f"RNA{RNA}_annotation", "interactions"]]
        small_df = small_df.groupby([f"RNA{RNA}_annotation"])
        if chimeras_or_fragments == "chimeras":
            small_df = small_df.count()
        elif chimeras_or_fragments == "fragments":
            small_df = small_df.sum()
        indexes = small_df.index
        for index in indexes:
            amount_dic[index] += small_df["interactions"][index]
    amount_list = [amount_dic[i] for i in options]
    amount_sum = sum(amount_list)
    return [i/amount_sum for i in amount_list]


def chimeras_annotations(chr_dic, experiments, base_path, annotation_path, rna_types_excel, chimeras_or_fragments):
    """
    Creates a bar plot of the annotation fraction for within and between chromosomes for each time point.
    :param chimeras_or_fragments: if "chimeras" create the bar plot according to the amount of chimeras. Else, according
    to the amount of fragments.
    """
    genes = get_annotation(annotation_path, separate_id_name=True)["name"].values.tolist()
    sRNAs, tRNAs, oRNAs = get_RNAs_types(rna_types_excel)
    options = ["3UTR", "sRNA", "CDS", "5UTR", "IGR", "AS", "tRNA", "oRNA", "IGT"]
    final_dic = {}
    chr1, chr2 = chr_dic.keys()
    chr1_name = chr_dic[chr1]
    chr2_name = chr_dic[chr2]
    for experiment in experiments:
        df = pd.read_excel(os.path.join(base_path, "RILseq_unified_results.xlsx"), sheet_name=experiment)
        df["RNA1_annotation"] = df["RNA1 name"].apply(get_gene_annotation, genes=genes, tRNAs=tRNAs, sRNAs=sRNAs, oRNAs=oRNAs)
        df["RNA2_annotation"] = df["RNA2 name"].apply(get_gene_annotation, genes=genes, tRNAs=tRNAs, sRNAs=sRNAs, oRNAs=oRNAs)
        chr1_df = df[(df["RNA1 chromosome"] == chr1) & (df["RNA2 chromosome"] == chr1)]
        chr2_df = df[(df["RNA1 chromosome"] == chr2) & (df["RNA2 chromosome"] == chr2)]
        chr1_chr2 = df[((df["RNA1 chromosome"] == chr1) & (df["RNA2 chromosome"] == chr2)) |
                           ((df["RNA1 chromosome"] == chr2) & (df["RNA2 chromosome"] == chr1))]
        if chimeras_or_fragments == "chimeras":
            chr1_amount = chr1_df.shape[0]
            chr1_chr2_amount = chr1_chr2.shape[0]
            chr2_amount = chr2_df.shape[0]
        else:
            chr1_amount = chr1_df["interactions"].sum()
            chr1_chr2_amount = chr1_chr2["interactions"].sum()
            chr2_amount = chr2_df["interactions"].sum()
        final_dic[f"{chr1_name}-{chr1_name} {experiment} ({chr1_amount})"] = get_RNAs_type_fractions(chr1_df, options, chimeras_or_fragments)
        final_dic[f"{chr1_name}-{chr2_name} {experiment} ({chr1_chr2_amount})"] = get_RNAs_type_fractions(chr1_chr2, options, chimeras_or_fragments)
        final_dic[f"{chr2_name}-{chr2_name} {experiment} ({chr2_amount})"] = get_RNAs_type_fractions(chr2_df, options, chimeras_or_fragments)
    df = pd.DataFrame.from_dict(final_dic, orient="index", columns=options)
    df = df[df.columns[(df != 0).any()]]
    fig = plt.figure()
    h = df.plot(kind="bar", linewidth=0, stacked=True, grid=False, cmap="tab20c")
    plt.subplots_adjust(bottom=0.35, right=0.7)
    plt.xticks(rotation=-45, ha='left')
    handles, labels = h.get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1], loc=(1.1, 0))
    plt.savefig(os.path.join(base_path, f"{chimeras_or_fragments}_annotation_fractions.png"))


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("config", help="Path to yaml configuration file.")

    args = parser.parse_args()
    config = load_config(args.config)

    for type_ in ("fragments", "chimeras"):
        for chrom in config["chr_dic"].keys():
            RNA1_RNA2_annotations(config["chr_dic"], config["base_path"], config["experiments"], config["rna_types_excel"], config["annotation_path"], type_, chrom)
        RNA1_RNA2_annotations(config["chr_dic"], config["base_path"], config["experiments"], config["rna_types_excel"], config["annotation_path"], type_)

        chimeras_annotations(config["chr_dic"], config["experiments"], config["base_path"], config["annotation_path"], config["rna_types_excel"], type_)


# if __name__ == '__main__':
#     main()




