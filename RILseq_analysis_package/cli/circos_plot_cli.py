import argparse
from RILseq_analysis_package import circos_plot
from RILseq_analysis_package.utils import load_config
import os

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "plot_type",
        choices=["two_genomes_half_circle", "two_genomes_real_proportions"],
        help="Which circos plot to generate"
    )
    parser.add_argument("config")
    parser.add_argument("--mark_step1", default=200000, help="interval of the scale marks of the first chromosome")
    parser.add_argument("--mark_step2", default=48500, help="interval of the scale marks of the second chromosome")
    parser.add_argument("--gene_list", default=None,
                        help="Plot only the chimeras of the genes in the list.")
    parser.add_argument("--present_chr", default=None, help="Plot only the chimeras of the given chromosome")
    args = parser.parse_args()
    config = load_config(args.config)

    base_path = config["base_path"]
    circos_path = os.path.join(base_path, "circos_plots")
    if not os.path.exists(circos_path):
        os.makedirs(circos_path)

    if args.plot_type == "two_genomes_half_circle":
        circos_plot.two_genomes_circos_plot_half_circle(circos_path, config["chr_sizes"], config["chr_dic"], args.present_chr, args.gene_list)
    elif args.plot_type == "two_genomes_real_proportions":
        circos_plot.two_genomes_circos_plot_real_proportions(circos_path, config["chr_sizes"], config["chr_dic"], args.mark_step1, args.mark_step2, args.present_chr)
