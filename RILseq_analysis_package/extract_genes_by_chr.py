import pandas as pd
from RILseq_analysis_package.utils import get_annotation, load_config
import argparse


def filter_chr_genes(annotation_path, old_file_path, new_file, chromosome, sep="\t", names_col="Unnamed: 0", identifier="name"):
    if identifier in ("name", "id"):
        annotations_df = get_annotation(annotation_path, chromosome=chromosome, separate_id_name=True)
        chr_genes = annotations_df[identifier].values.tolist()
    else:
        annotations_df = get_annotation(annotation_path, chromosome=chromosome, identifier=identifier)
        chr_genes = annotations_df["identifier"].values.tolist()

    old_file_df = pd.read_csv(old_file_path, sep=sep)

    new_df = old_file_df[old_file_df[names_col].isin(chr_genes)]
    new_df.to_csv(new_file, index=False, sep=sep)


def main():
    parser = argparse.ArgumentParser(description="Filter a file and keep only genes belonging to a specified chromosome.")
    parser.add_argument("input", help="Path to the input file.")
    parser.add_argument("output", help="Path to the output file.")
    parser.add_argument("chr", help="Chromosome to keep, according to the chromosome names in the annotation file.")
    parser.add_argument("config", help="Path to the YAML configuration file.")
    parser.add_argument("--sep", default="\t", help="Column separator used in the input and output files (default: tab).")
    parser.add_argument("--names_col", default="Unnamed: 0", help="Name of the column containing gene names (default: there is no column name).")
    parser.add_argument("--identifier", default="name", help="Type of gene identifier to use from the annotation file (default: name).")

    args = parser.parse_args()

    filter_chr_genes(load_config(args.config)["annotation_path"], args.input, args.output, args.chr, args.sep, args.names_col, args.identifier)
