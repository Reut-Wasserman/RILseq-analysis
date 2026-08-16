import argparse
from RILseq_analysis_package.generate_and_edit_RILseq_xslx import add_number_of_libraries
from RILseq_analysis_package.utils import load_config


def main():
    parser = argparse.ArgumentParser(description="Record for each chimera the number of libraries in which it appears, "
                                                 "and rename columns in the S-chimeras workbook.")

    parser.add_argument("config", help="Path to the YAML configuration file.")

    args = parser.parse_args()
    config = load_config(args.config)
    add_number_of_libraries(config["base_path"], config["experiments"], config["replicates"], config["chr_dic"])

