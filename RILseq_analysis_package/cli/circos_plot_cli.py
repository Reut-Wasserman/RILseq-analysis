import argparse
from RILseq_analysis_package import circos_plot
from RILseq_analysis_package.utils import load_config
import os

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("plot_type", choices=["two_genomes_half_circle", "two_genomes_real_proportions"], help="Which circos plot to generate")
    parser.add_argument("config")
    parser.add_argument("--mark_step1", default=200000, help="interval of the scale marks of the first chromosome")
    parser.add_argument("--mark_step2", default=48500, help="interval of the scale marks of the second chromosome")
    parser.add_argument("--genes", default=None,
                        help="Plot only the chimeras of the genes in the list.")
    parser.add_argument("--present_chr", default=None, help="Plot only the chimeras of the given chromosome")
    parser.add_argument("--out_path", default=None)
    parser.add_argument("--experiments", default=None, help="list of experiments to create the files for them, seperated with a comma. "
                                                            "The experiments must be as the sheet names in the RILseq_unified_results.excel file")
    args = parser.parse_args()
    config = load_config(args.config)

    base_path = config["base_path"]
    if args.out_path is None:
        circos_path = os.path.join(base_path, "circos_plots")
    else:
        circos_path = args.out_path

    if not os.path.exists(circos_path):
        os.makedirs(circos_path)

    if args.genes is not None:
        genes_list = args.genes.split("_")
    else:
        genes_list = None

    if args.experiments is not None:
        experiments = args.experiments.split(",")
    else:
        experiments = None

    if args.plot_type == "two_genomes_half_circle":
        circos_plot.two_genomes_circos_plot_half_circle(circos_path, base_path, config["chr_len"], config["chr_dic"], args.present_chr, genes_list, experiments)
    elif args.plot_type == "two_genomes_real_proportions":
        circos_plot.two_genomes_circos_plot_real_proportions(circos_path, base_path, config["chr_len"], config["chr_dic"], args.mark_step1, args.mark_step2, args.present_chr, genes_list)


