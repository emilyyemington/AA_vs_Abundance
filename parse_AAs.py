import re
from collections import Counter
import pandas as pd


def split_ASTK(c1, c2, c3, astk):
    frs = re.split(c1 + '|' + c2 + '|' + c3, astk)
    return frs


def count_AAs(sequence, reads) -> Counter:
    aa_counter = Counter()
    for AA in sequence:
        aa_counter[AA] += float(reads)
    return aa_counter


def parse_AAs(summ_data, sequence_data):
    summary_data = summ_data.copy()
    seq_data = sequence_data.copy()

    # get FRs from CDRs and ASTK

    FR_data = seq_data.apply(lambda x:
                             split_ASTK(x['CDRH1_aa'], x['CDRH2_aa'], x['CDRH3_aa'], x['ASTK']),
                             axis=1,
                             result_type='expand')
    FR_data.columns = ["FR1", "FR2", "FR3", "FR4"]

    seq_data = pd.concat([seq_data, FR_data], axis=1)

    # Time to count the amino acids
    AA_data = pd.DataFrame(seq_data["ClusterID"].copy())
    AA_data.columns = ["Clusters"]

    # counts for every sequence type and adds to new columns
    seq_types = ["CDRH1_aa", "CDRH2_aa", "CDRH3_aa", "FR1", "FR2", "FR3", "FR4"]
    for seq in seq_types:
        AA_data[seq + "_count"] = pd.DataFrame(seq_data.apply(lambda x: count_AAs(x[seq], x['Collapsed']), axis=1))
    # sums up each sequence type for overall total counts
    AA_data["total_CDR_count"] = AA_data["CDRH1_aa_count"] + AA_data["CDRH2_aa_count"] + AA_data["CDRH3_aa_count"]
    AA_data["total_FR_count"] = AA_data["FR1_count"] + AA_data["FR2_count"] + AA_data["FR3_count"]
    AA_data["total_count"] = AA_data["total_CDR_count"] + AA_data["total_FR_count"]

    # sums up counts by cluster ID
    AA_data = AA_data.groupby("Clusters").sum().copy()
    summary_data = summary_data.join(AA_data, on="Clusters")

    seq_types = seq_types + ["total_CDR", "total_FR", "total"]  # for column names

    # normalizes AA counts
    for seq in seq_types:
        summary_data[seq + "_norm"] = summary_data.apply(
            lambda x: {key: round(x[seq + "_count"][key] / x["VH Reads"], 2) for key in x[seq + "_count"]}, axis=1)

    return summary_data, seq_types
