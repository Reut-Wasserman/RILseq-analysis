import argparse
from RILseq_analysis_package import circos_plot
from RILseq_analysis_package.utils import load_config
import os

def main():
    parser = argparse.ArgumentParser(description="Create a Circos plot showing RIL-seq chimeric interactions between two chromosomes.")

    parser.add_argument("plot_type", choices=["two_genomes_half_circle", "two_genomes_real_proportions"], help="Type of circos plot to generate.")
    parser.add_argument("config", help="Path to the YAML configuration file.")
    parser.add_argument("--mark_step1", default=200000, type=int, help="Distance between scale marks on the first chromosome (default: 200000).")
    parser.add_argument("--mark_step2", default=48500, type=int, help="Distance between scale marks on the second chromosome (default: 48500).")
    parser.add_argument("--genes", default=None,
                        help="Plot only chimeras involving the specified genes. Separate gene names with underscores.")
    parser.add_argument("--present_chr", default=None, help="Plot only chimeras involving the specified chromosome.")
    parser.add_argument("--out_path", default=None,
                        help="Output directory for the Circos plot. If not provided, a 'circos_plots' directory is created under base_path.")
    parser.add_argument("--experiments", default=None,
                        help="Comma-separated list of experiment names to plot. "
                             "Experiment names must match sheet names in RILseq_unified_results.xlsx.")
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


