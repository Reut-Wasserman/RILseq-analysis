import argparse
from RILseq_analysis_package.generate_and_edit_RILseq_xslx import merge_RILseq_results
from RILseq_analysis_package.utils import load_config


def main():
    parser = argparse.ArgumentParser(description="Merge RIL-seq results from multiple experiments into a unified Excel file.")

    parser.add_argument("config", help="Path to the YAML configuration file.")
    parser.add_argument("file_names",
                        help="path to file contains two columns. One is the name of the sig_interactions file names and "
                             "the other is the corresponding experiment names as appear in the yaml file. The two columns "
                             "are separated by a Tab.")
    parser.add_argument("--no-annotations", dest="annotations", action="store_false", default=True,
                        help="Do not add RNA annotations to the merged results.")

    args = parser.parse_args()
    config = load_config(args.config)
    merge_RILseq_results(config["base_path"], config["annotation_path"], config["rna_types_excel"], args.file_names, args.annotations)

