import os
import pandas as pd
from RILseq_analysis_package.utils import get_annotation


def convert_to_kb_mb(range_lst):
    if range_lst[-1] > 10 ** 6:
        return [str("%.0f" % (i / 10 ** 6)) + " MB" for i in range_lst]
    else:
        return [str("%.0f" % (i / 1000)) + " KB" for i in range_lst]


def get_scale_lists(chr_sizes, chr_dic, factor, kb_mb):
    """
    Returns lists for the scale marks of the circos plot, the multiplication factor needed for the shorter chromosome, and its name.
    """
    # chr_sizes = {}
    # with open(CHR_SIZES_PATH, "r") as f:
    #     lines = f.readlines()
    #     max_len = 0
    #     min_len = 10**9
    #     min_chr = ""
    #     for line in lines:
    #         chrom, chr_len = line.split()
    #         chr_sizes[chrom] = chr_len
    #         if int(chr_len) > max_len:
    #             max_len = int(chr_len)
    #         if int(chr_len) < min_len:
    #             min_len = int(chr_len)
    #             min_chr = chrom
    max_chr = max(chr_sizes, key=chr_sizes.get)
    min_chr = min(chr_sizes, key=chr_sizes.get)

    max_len = chr_sizes[max_chr]
    min_len = chr_sizes[min_chr]

    multiply_factor = round(max_len/min_len)

    ranges = {}
    for chrom in chr_sizes.keys():
        ranges[chrom] = list(range(0, int(chr_sizes[chrom]), 10 ** (len(chr_sizes[chrom]) - factor)))

    chromosome = []
    for i in ranges.keys():
        chromosome += [chr_dic[i]] * len(ranges[i])

    chrom_start = []
    for i in ranges.keys():
        if i == min_chr:
            chrom_start += [i * multiply_factor for i in ranges[i]]
        else:
            chrom_start += ranges[i]

    name = []
    if kb_mb:
        for chrom in chr_sizes.keys():
            name += convert_to_kb_mb(ranges[chrom])
    else:
        for i in ranges.keys():
            name += [str(k) for k in ranges[i]]

    return chromosome, chrom_start, name, multiply_factor, min_chr


def add_color(chimeras_df, chr1, chr2):
    """
    Add a column of colors to the chimeras. Set different colors to chimeras from different genomes and to chimeras between genomes.
    """
    chimeras_df["PlotColor"] = "."
    chimeras_df["PlotColor"][(chimeras_df["RNA1 chromosome"] == chr1) & (chimeras_df["RNA2 chromosome"] == chr1)] = "gray48"
    chimeras_df["PlotColor"][(chimeras_df["RNA1 chromosome"] == chr1) & (chimeras_df["RNA2 chromosome"] == chr2)] = "blue"
    chimeras_df["PlotColor"][(chimeras_df["RNA1 chromosome"] == chr2) & (chimeras_df["RNA2 chromosome"] == chr1)] = "blue"
    chimeras_df["PlotColor"][(chimeras_df["RNA1 chromosome"] == chr2) & (chimeras_df["RNA2 chromosome"] == chr2)] = "red"
    return chimeras_df


def two_genomes_circos_plot_half_circle(base_path, chr_sizes, chr_dic, present_chr=None, gene_list=None):
    """
    Generates files required for creating circos plot of two genomes, each represented by half a circle.
     Use the output files for circos_plot.R script.
    :param gene_list: list of genes. If it's not None the circos plot will contain only the chimeras of the genes in the list.
    """
    chromosome, chrom_start, gene, multiply_factor, chr_to_multiply = get_scale_lists(chr_sizes, chr_dic, 2, kb_mb=False)
    scale_df = pd.DataFrame({"Chromosome":chromosome, "chromStart":chrom_start, "chromEnd":[i+1 for i in chrom_start],
                             "Gene":gene})

    scale_df.to_csv(os.path.join(base_path, "circos_plots", "two_genomes_scale_each_half_circle.csv"), index=False)

    chromosome, chrom_start, gene, _, _ = get_scale_lists(chr_sizes, chr_dic, 1, kb_mb=True)

    scale_df = pd.DataFrame({"Chromosome":chromosome, "chromStart":chrom_start, "chromEnd":[i+1 for i in chrom_start],
                             "Gene":gene})
    scale_df.to_csv(os.path.join(base_path, "circos_plots", "two_genomes_scale_numbers_each_half_circle.csv"), index=False)

    RILseq_excel = pd.ExcelFile(os.path.join(base_path, "RILseq_unified_results.xlsx"))
    for experiment in RILseq_excel.sheet_names:
        chimeras_df = RILseq_excel.parse(experiment)
        if gene_list is not None:
            chimeras_df = chimeras_df[(chimeras_df["RNA1 name"].isin(gene_list)) | (chimeras_df["RNA2 name"].isin(gene_list))]
        if present_chr:
            chimeras_df = chimeras_df[(chimeras_df["RNA1 chromosome"] == present_chr) & (chimeras_df["RNA2 chromosome"] == present_chr)]

        for i in ("1", "2"):
            for j in [f"Start of RNA{i} first read", f"Start of RNA{i} last read"]:
                chimeras_df[j][chimeras_df[f"RNA{i} chromosome"] == chr_to_multiply] *= multiply_factor
        chr1, chr2 = chr_sizes.keys
        chimeras_df = add_color(chimeras_df, chr1, chr2)

        for genome, genome_name in chr_dic.items():
            chimeras_df["RNA1 chromosome"][chimeras_df["RNA1 chromosome"] == genome] = genome_name
            chimeras_df["RNA2 chromosome"][chimeras_df["RNA2 chromosome"] == genome] = genome_name

        file_name = f"{experiment}_two_genomes_chimeras_each_half_circle"
        if gene_list is not None:
            genes = "_".join(gene_list)
            file_name = f"{genes}_chimeras_{experiment}"
        chimeras_df.to_csv(os.path.join(base_path, "circos_plots", f"{file_name}.csv"),
                           columns=["RNA1 chromosome", "Start of RNA1 first read", "Start of RNA1 last read", "RNA2 chromosome", "Start of RNA2 last read", "Start of RNA2 first read", "PlotColor"],
                           header=["Chromosome", "chromStart", "chromEnd", "Chromosome.1", "chromStart.1", "chromEnd.1", "PlotColor"], index=False)


def two_genomes_circos_plot_real_proportions(base_path, chr_sizes, chr_dic, mark_step1, mark_step2, present_chr=None):
    """
    Generates files required for creating circos plot of two genomes. Use the output files for circos_plot.R script.
    """
    # annotations_df = get_annotation(separate_id_name=True)
    # annotations_df.rename(columns={"name": "Gene", "chr": "Chromosome", "start": "chromStart", "end": "chromEnd"}, inplace=True)
    # for genome, genome_name in CHR_DIC.items():
    #     annotations_df["Chromosome"][annotations_df["Chromosome"] == genome] = genome_name
    # annotations_df.to_csv(os.path.join(BASE_PATH, r"RILseq\circos_plots\two_genomes_genes_real_proportions.csv"),
    #                        columns=["Chromosome", "chromStart", "chromEnd", "Gene"],
    #                        header=["Chromosome", "chromStart", "chromEnd", "Gene"], index=False)

    chr1, chr2 = chr_sizes.keys
    chrI_range = list(range(0, chr_sizes[chr1], mark_step1))
    chrII_range = list(range(0, chr_sizes[chr2], mark_step2))
    scale_df = pd.DataFrame({"Chromosome":[chr_sizes[chr1]]*len(chrI_range) + [chr_sizes[chr2]]*len(chrII_range),
                             "chromStart":chrI_range + chrII_range,
                             "chromEnd":[i+1 for i in chrI_range] +[i+1 for i in chrII_range],
                             "Gene":convert_to_kb_mb(chrI_range+chrII_range)})
    scale_df.to_csv(os.path.join(base_path, "circos_plots", "two_genomes_scale_real_proportions.csv"), index=False)

    RILseq_excel = pd.ExcelFile(os.path.join(base_path, "RILseq_unified_results.xlsx"))
    for experiment in RILseq_excel.sheet_names:
        chimeras_df = RILseq_excel.parse(experiment)

        if present_chr:
            chimeras_df = chimeras_df[(chimeras_df["RNA1 chromosome"] == present_chr) & (chimeras_df["RNA2 chromosome"] == present_chr)]

        chimeras_df = add_color(chimeras_df, chr1, chr2)

        for genome, genome_name in chr_dic.items():
            chimeras_df["RNA1 chromosome"][chimeras_df["RNA1 chromosome"] == genome] = genome_name
            chimeras_df["RNA2 chromosome"][chimeras_df["RNA2 chromosome"] == genome] = genome_name

        chimeras_df.to_csv(os.path.join(base_path, "circos_plots", f"{experiment}_two_genomes_chimeras_real_proportions.csv"),
                           columns=["RNA1 chromosome", "Start of RNA1 first read", "Start of RNA1 last read", "RNA2 chromosome", "Start of RNA2 last read", "Start of RNA2 first read", "PlotColor"],
                           header=["Chromosome", "chromStart", "chromEnd", "Chromosome.1", "chromStart.1", "chromEnd.1", "PlotColor"], index=False)



# if __name__ == "__main__":
#     two_genomes_circos_plot_half_circle()
#     two_genomes_circos_plot_half_circle(["preS"])
#     two_genomes_circos_plot_half_circle(["lpr2"])
#     two_genomes_circos_plot_real_proportions()
#     E_coli_circos_plot()



