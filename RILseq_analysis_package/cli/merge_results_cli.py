import argparse
from RILseq_analysis_package.generate_and_edit_RILseq_xslx import merge_RILseq_results
from RILseq_analysis_package.utils import load_config


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("config")
    parser.add_argument("--annotations", default=True)

    args = parser.parse_args()
    config = load_config(args.config)
    merge_RILseq_results(config["base_path"], config["annotation_path"], config["rna_types_excel"], args.annotations)

