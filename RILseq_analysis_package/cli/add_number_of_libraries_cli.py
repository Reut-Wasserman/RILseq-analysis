import argparse
from RILseq_analysis_package.generate_and_edit_RILseq_xslx import add_number_of_libraries
from RILseq_analysis_package.utils import load_config


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("config")

    args = parser.parse_args()
    config = load_config(args.config)
    add_number_of_libraries(config["base_path"], config["experiments"], config["replicates"], config["chr_dic"])

