import os
import pandas as pd
from RILseq_analysis_package.utils import get_annotation, get_RNA_types, get_replicates
# from RILseq_analysis_package.defaults import *


def find_genomic_annotation(gene, sRNAs_list, tRNAs_list, all_genes_list):
    if gene in sRNAs_list:
        return "sRNA"
    if gene in tRNAs_list:
        return "tRNA"
    if "3UTR" in gene:
        return "3UTR"
    if "5UTR" in gene:
        return "5UTR"
    if "AS" in gene:
        return "AS"
    if "IGR" in gene:
        return "IGR"
    if "IGT" in gene:
        return "IGT"
    if gene in all_genes_list:
        return "CDS"
    return "unknown"


def merge_RILseq_results(base_path, annotation_path, rna_types_excel, add_genomic_annotation=True):
    unified = pd.ExcelWriter(os.path.join(base_path, 'RILseq_unified_results.xlsx'))
    single = pd.ExcelWriter(os.path.join(base_path, 'RILseq_single_results.xlsx'))
    for file in os.listdir(base_path):
        if file.endswith("sig_interactions.txt"):
            df = pd.read_csv(os.path.join(base_path, file), sep="\t")
            if not df.empty:
                if add_genomic_annotation:
                    sRNAs_list = get_RNA_types("sRNA", rna_types_excel)
                    tRNAs_list = get_RNA_types("tRNA", rna_types_excel)
                    genes_names = get_annotation(annotation_path, separate_id_name=True)["name"].values.tolist()
                    df["Genomic annotation of RNA1"] = df["RNA1 name"].apply(find_genomic_annotation, sRNAs_list=sRNAs_list, tRNAs_list=tRNAs_list, all_genes_list=genes_names)
                    df["Genomic annotation of RNA2"] = df["RNA2 name"].apply(find_genomic_annotation, sRNAs_list=sRNAs_list, tRNAs_list=tRNAs_list, all_genes_list=genes_names)
                if file.startswith("unified"):
                    sheet_name = file.replace("unified_", "").replace("_all_fragments_l25.txt_sig_interactions.txt", "").replace("_mapping", "")
                    df.to_excel(unified, sheet_name, index=False)
                else:
                    sheet_name = file.replace("cutadapt_bwa.bam_mapping_all_fragments_l25.txt_sig_interactions.txt", "S_chimeras")
                    sheet_name = sheet_name.replace("cutadapt_bwa.bam_sig_interactions.txt", "S_chimeras")
                    if sheet_name.startswith("RILSeq_"):
                        sheet_name = sheet_name.replace("RILSeq_", "")
                    df.to_excel(single, sheet_name, index=False)
    unified.save()
    single.save()


def find_number_of_libraries_helper(name1, name2, start1, end1, start2, end2, chr1, chr2, strand1, strand2, df):
    cur1 = df[(df["RNA1 name"] == name1) & (df["RNA2 name"] == name2) & (df["RNA1 strand"] == strand1) & (df["RNA2 strand"] == strand2) & (df["RNA1 chromosome"] == chr1) & (df["RNA2 chromosome"] == chr2)]
    cur1 = cur1[~((cur1["Start of RNA1 first read"] > end1) | (cur1["Start of RNA1 last read"] < start1))]
    cur1 = cur1[~((cur1["Start of RNA2 last read"] > end2) | (cur1["Start of RNA2 first read"] < start2))]
    return cur1.shape[0]


def find_number_of_libraries(unify_chimera, singles):
    name1, name2, start1, end1, start2, end2, chr1, chr2, strand1, strand2 = unify_chimera[["RNA1 name", "RNA2 name", "Start of RNA1 first read", "Start of RNA1 last read", "Start of RNA2 last read", "Start of RNA2 first read", "RNA1 chromosome", "RNA2 chromosome", "RNA1 strand", "RNA2 strand"]]

    counter = 0
    for single in singles:
        single = single.astype({"Start of RNA1 first read":int, "Start of RNA2 first read":int, "Start of RNA1 last read":int, "Start of RNA2 last read":int})
        amount1 = find_number_of_libraries_helper(name1, name2, start1, end1, start2, end2, chr1, chr2, strand1, strand2, single)
        amount2 = find_number_of_libraries_helper(name2, name1, start2, end2, start1, end1, chr2, chr1, strand2, strand1, single)
        if amount1 + amount2 != 0:
            counter += 1

    if counter == 0:
        return "U"
    return counter


def add_number_of_libraries(base_path, experiments, replicates, chr_dic):
    new = pd.ExcelWriter(os.path.join(base_path, 'RILseq_unified_results_with_number_of_libraries.xlsx'))
    single_results_excel = pd.ExcelFile(os.path.join(base_path, "RILseq_single_results.xlsx"))
    for experiment in experiments:
        singles = []
        for i in replicates:
            replicate_name = f"{experiment + i}_S_chimeras"
            if replicate_name in single_results_excel.sheet_names:
                single1 = single_results_excel.parse(f"{replicate_name}_S_chimeras")
                singles.append(single1)

        unify_df = pd.read_excel(os.path.join(base_path, "RILseq_unified_results.xlsx"), sheet_name=experiment)

        unify_df["# of libraries"] = unify_df.apply(find_number_of_libraries, singles=singles, axis=1)
        unify_df = unify_df[["RNA1 name", "RNA2 name", "interactions", "# of libraries", "Normalized Odds Ratio (NOR)", "odds ratio", "Fisher's exact test p-value", "Genomic annotation of RNA1", "Genomic annotation of RNA2", "RNA1 description", "RNA2 description", "RNA1 chromosome", "Start of RNA1 first read", "Start of RNA1 last read", "RNA1 strand", "RNA2 chromosome", "Start of RNA2 last read", "Start of RNA2 first read", "RNA2 strand", "other interactions of RNA1", "other interactions of RNA2", "total other interactions", "total RNA reads1", "total RNA reads2", "lib norm IP RNA1", "lib norm IP RNA2", "lib norm total RNA1", "lib norm total RNA2", "IP/total ratio1", "IP/total ratio2", "RNA1 EcoCyc ID", "RNA2 EcoCyc ID"]]
        unify_df.rename(columns={"interactions":"# of chimeric fragments", "Normalized Odds Ratio (NOR)":"Normalized Odds Ratio", "odds ratio":"Odds Ratio",
                         "Start of RNA1 first read":"RNA1 from", "Start of RNA1 last read":"RNA1 to", "Start of RNA2 last read":"RNA2 from",
                         "Start of RNA2 first read":"RNA2 to", "other interactions of RNA1":"other fragments of RNA1", "other interactions of RNA2":"other fragments of RNA2",
                         "total other interactions":"Total other fragments", "total RNA reads1":"RNA1 in total RNA (# of reads)",
                         "total RNA reads2":"RNA2 in total RNA (# of reads)"}, inplace=True)
        unify_df["RNA1 chromosome"] = unify_df["RNA1 chromosome"].apply(lambda x: chr_dic[x])
        unify_df["RNA2 chromosome"] = unify_df["RNA2 chromosome"].apply(lambda x: chr_dic[x])
        unify_df.to_excel(new, experiment, index=False)

    new.save()


# if __name__ == '__main__':
#     merge_RILseq_results(rf"{BASE_PATH}\RILSeq\results", add_genomic_annotation=True)
#     add_number_of_libraries(rf"{BASE_PATH}\RILSeq\results", get_experiments())