import pandas as pd
import yaml
# from RILseq_analysis_package.defaults import *


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def get_gene_identifier(gene_id, identifier):
    """
    For the given gene_id, returns the value of the given identifier.
    """
    split_id = gene_id.split(";")
    if identifier in gene_id:
        for i in split_id:
            if identifier in i:
                if i.startswith(" "):
                    return i.split(" ")[2].replace("\"", "")
                return i.split(" ")[1].replace("\"", "")
    id_ = split_id[0].split(" ")
    return id_[1].replace("\"", "")


def get_annotation(annotation_path, chromosome=None, separate_id_name=False, identifier=None):
    """
    Returns a df of the annotations according to ANNOTATION_FILE.
    separate_id_name: if True add a column of name, the id column will contain only the id. Use when there are two
                      identifiers - id and name.
    identifier: add a column "identifier" which contains the values of the give identifier. Use when there are other
                identifier (not id and name).
    """
    cols_names = ["chr", "EcoCyc", "exon", "start", "end", ".1", "strand", ".2", "id"]
    df = pd.read_csv(annotation_path, sep="\t", names=cols_names)
    if chromosome:
        if type(chromosome) is str:
            df = df[df["chr"] == chromosome]
        if type(chromosome) is list:
            df = df[df["chr"].isin(chromosome)]
    if separate_id_name:
        df["name"] = df["id"].apply(lambda x: x.split()[3].replace(";", "").replace("\"", ""))
        df["id"] = df["id"].apply(lambda x: x.split()[1].replace(";", "").replace("\"", ""))
    if identifier:
        df["identifier"] = df["id"].apply(get_gene_identifier, args=(identifier,))
    return df


def get_only_id(gene):
    return gene.split()[1].replace(";", "").replace("\"", "")


def get_only_name(gene):
    return gene.split()[3].replace(";", "").replace("\"", "")


def get_sRNAs(rna_types_excel, organism=None):
    if organism == "lambda":
        return ["preS", "lpr2", "6S RNA"]
    sRNAs = pd.read_excel(rna_types_excel, sheet_name="sRNA")
    sRNAs = sRNAs["Name"].values.tolist()
    if EXPERIMENT == "lambda" and organism is None:
        return sRNAs + ["preS", "lpr2", "6S RNA"]
    return sRNAs


def get_RNA_types(RNA_type, rna_types_excel):
    if RNA_type == "sRNA":
        return get_sRNAs(rna_types_excel)
    RNAs = pd.read_excel(rna_types_excel, sheet_name=RNA_type)
    RNAs = RNAs["Name"].values.tolist()
    if RNA_type == "oRNA":
        RNAs += ["pspH", "RirA"]
    return RNAs


def remove_UTR(df):
    df["RNA1 name"] = df["RNA1 name"].apply(lambda x: x.split(".EST3UTR")[0].split(".EST5UTR")[0].split(".5UTR")[0].split(".3UTR")[0])
    df["RNA2 name"] = df["RNA2 name"].apply(lambda x: x.split(".EST3UTR")[0].split(".EST5UTR")[0].split(".5UTR")[0].split(".3UTR")[0])
    return df


def get_gene_chr(gene, annotation_df):
    return annotation_df["chr"][annotation_df["name"] == gene]


def get_E_coli_lambda_experiments(with_wt):
    res = []
    for time in ("30", "60"):
        for experiment in ("Hfq_", "Hfq_lambda_", "wt_lambda_", "wt"):
            if not with_wt and "wt" in experiment:
                continue
            if "wt" in experiment:
                if time == "60":
                    continue
            res.append(experiment + time)
    return res


def get_experiments(with_wt=True):
    if EXPERIMENT == "lambda":
        return get_E_coli_lambda_experiments(with_wt)


def get_replicates():
    if EXPERIMENT == "lambda":
        return ["_I", "_II"]

print(get_E_coli_lambda_experiments(True))
