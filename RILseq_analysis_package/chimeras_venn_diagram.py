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
              set_colors=("SkyBlue", "Salmon"))
        plt.savefig(os.path.join(output_path, f"{conditions[0]}_and_{conditions[1]}_venn.jpg"))
        with open(os.path.join(output_path, f"{conditions[0]}_and_{conditions[1]}_chimeras_pairs.txt"), "w") as f:
            pairs = experiments_chimeras[conditions[0]] & experiments_chimeras[conditions[1]]
            for i in pairs:
                f.write(i[0] + "," + i[1] + "\n")


def main():
    parser = argparse.ArgumentParser("Generate Venn diagrams of chimeras pairs between two conditions. RNA1/RNA2 positions are ignored.")

    parser.add_argument("config", help="Path to the YAML configuration file.")
    parser.add_argument("conditions_pairs", help="Pairs of experiments to compare, with the two conditions separated by '-'. "
     "Multiple pairs can be separated by commas. "
     "Example: exp1-exp2,exp3-exp4 creates two Venn diagrams: one comparing exp1 and exp2, and the other comparing exp3 and exp4.")
    parser.add_argument("--chr", help="Comma-separated list of chromosomes to include. If not provided, chimeras from all chromosomes are included.")

    args = parser.parse_args()

    base_path = load_config(args.config)["base_path"]
    venn_path = os.path.join(base_path, "venn_diagrams")
    if not os.path.exists(venn_path):
        os.makedirs(venn_path)

    pairs = args.conditions_pairs.split(",")
    pairs = [i.split("-") for i in pairs]

    if args.chr is not None:
        chr_list = args.chr.split(",")
    else:
        chr_list = None
    plot(chr_list, pairs, venn_path, base_path)

