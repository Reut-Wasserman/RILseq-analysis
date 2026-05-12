import argparse
import os
import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from RILseq_analysis_package.utils import load_config


def get_chimeras_pairs(RILseq_path, experiment, chromosomes=None):
    """
    Reads the RILSeq results and returns a set of the chimeras pairs. Doesn't take into account which RNA is first and
    which is second in the chimera.
    """
    df = pd.read_excel(os.path.join(RILseq_path, "RILseq_unified_results.xlsx"), sheet_name=experiment)
    if chromosomes:
        df = df[(df["RNA1 chromosome"].isin(chromosomes)) & (df["RNA2 chromosome"].isin(chromosomes))]
    chimeras_pairs = df[["RNA1 name", "RNA2 name"]].values
    chimeras_set = set()
    for pair in chimeras_pairs:
        rna1, rna2 = pair
        if rna1 > rna2:
            chimeras_set.add((rna1, rna2))
        else:
            chimeras_set.add((rna2, rna1))
    return chimeras_set


def get_all_conditions(conditions):
    all_conditions = []
    for conditions_tuple in conditions:
        all_conditions.append(conditions_tuple[0])
        all_conditions.append(conditions_tuple[1])
    return all_conditions


def plot(chr_lst, conditions_to_compare, output_path, RILseq_path):
    """Plot the venn diagrams according to VENN_CONDITIONS and writes the overlapping pairs into a file."""
    experiments_chimeras = {}
    for condition in get_all_conditions(conditions_to_compare):
        if chr_lst:
            chimeras_pairs = get_chimeras_pairs(RILseq_path, condition, chromosomes=chr_lst)
        else:
            chimeras_pairs = get_chimeras_pairs(RILseq_path, condition)
        experiments_chimeras.update({condition: chimeras_pairs})

    for conditions in conditions_to_compare:
        plt.figure()
        venn2([experiments_chimeras[conditions[0]], experiments_chimeras[conditions[1]]],
              set_labels=(conditions[0], conditions[1]),
              set_colors=("SkyBlue", "Salmon", "SkyBlue"))
        plt.savefig(os.path.join(output_path, f"{conditions[0]}_and_{conditions[1]}_venn.jpg"))
        with open(os.path.join(output_path, f"{conditions[0]}_and_{conditions[1]}_chimeras_pairs.txt"), "w") as f:
            pairs = experiments_chimeras[conditions[0]] & experiments_chimeras[conditions[1]]
            for i in pairs:
                f.write(i[0] + "," + i[1] + "\n")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("config", help="Path to yaml configuration file.")
    parser.add_argument("output_path", help="Output path.")
    parser.add_argument("conditions_pairs",
                        help="List of tuples. Each tuple contain two conditions to plot.")
    parser.add_argument("--chr", help="List of chromosomes. Include chimeras from those chromosomes only.")

    args = parser.parse_args()

    plot(args.chr, args.conditions_pairs, args.output_path, load_config(args.config)["base_path"])




