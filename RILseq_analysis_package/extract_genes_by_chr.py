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
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("chr")
    parser.add_argument("config")
    parser.add_argument("--sep", default="\t")
    parser.add_argument("--names_col", default="Unnamed: 0")
    parser.add_argument("--identifier", default="name")

    args = parser.parse_args()

    filter_chr_genes(load_config(args.config)["annotation_path"], args.input, args.output, args.chr, args.sep, args.names_col, args.identifier)


# if __name__ == '__main__':
#     file = os.path.join(BASE_PATH, r"RNASeq\all_counts_table.txt")
#     filter_chr_genes(file, os.path.join(BASE_PATH, r"RNASeq\all_counts_table_E_coli_genes.txt"), "chrI")
#     filter_chr_genes(file, os.path.join(BASE_PATH, r"RNASeq\all_counts_table_lambda_genes.txt"), "chrII")
