#
# This script takes in IgSeq files and counts residues of the CDRH and FR regions of the top ten most abundant mAbs,
# normalizing/averaging them by read count and outputting them into an Excel file by percent AA makeup out of the total
# length of the region. The outputted file will contain one summary sheet with all information, and several other
# sheets with the exact AA makeup of each region.
#
# For example, a mAb with a ten-residue-average CDRH1 region with an average of two serines will be S: 0.20 in a
# sheet called "CDRH1_aa_summary" in the outputted file
#
# The IgSeq files should be formatted such that:
#   -> they are Excel files
#   -> each file name contains the donor ID / day information
#   -> there are two sheets within each file, one named "Summary" and the other "VH Database"
#       -> in "Summary", there are columns labeled "Clusters", "Percent", "VH Reads", and "VH Gene"
#       -> "Summary" is either sorted by abundance OR only contains the top 10 clusters
#       -> in "VH Database" are columns for CDRH1, 2, and 3 sequences, labeled "CDRH#_aa" where # is 1, 2, or 3
#       -> in "VH Database" there is a full sequence containing the framework regions labeled "ASTK"
#       -> in "VH Database" there is a column with the number of reads labeled "Collapsed"
#       -> in "VH Database" there is a column labeled "ClusterID"

# For ease of choosing files, I would recommend putting them all in the same folder so that you can click on more than
# one at a time in the file selector

# region Import packages
import re
from collections import Counter
import pandas as pd
from tkinter import filedialog as fd
from parse_AAs import *

# endregion

# region File Import and Reading
files = fd.askopenfilenames(title="Select Formatted IgSeq Excel Files")
summary_data = pd.DataFrame()
seq_data = pd.DataFrame()
clusters = []

# Finds and saves data for the top ten clusters in every file
for file in files:
    # Reads relevant summary data for top 10 clusters from each file
    summ_cols = ['Clusters', 'Percent', 'VH Reads', 'VH Gene']
    file_data = pd.read_excel(file, sheet_name='Summary',
                              usecols=summ_cols).truncate(after=9)
    file_data["File"] = file  # keeps track of which donor has which cluster
    summary_data = pd.concat([summary_data, file_data], axis=0)

    # Reads sequencing data from the top 10 clusters
    clusters.extend(summary_data["Clusters"].tolist())
    seq_data = pd.read_excel(file, sheet_name='VH Database')

    # Checks if VH Database is formatted correctly
    if not all(x in seq_data.columns for x in ['CDRH1_aa', 'CDRH2_aa', 'CDRH3_aa', 'ASTK', 'ClusterID']):
        raise NameError("Check formatting of VH Database")

    seq_data = seq_data[seq_data['ClusterID'].isin(clusters)]  # only saves data for relevant clusters

# endregion

# region Write AA data to Excel

# Save summary to file name
filename = fd.asksaveasfilename(defaultextension=".xlsx")

with pd.ExcelWriter(filename) as writer:
    # get and reorganize data from AA parser
    summary_data, seq_types = parse_AAs(summary_data, seq_data)
    # save to summary sheet
    summary_data.to_excel(writer, sheet_name="Script Summary", index=True)

    # Save individual sheets for each sequence type
    seq_dict = {}
    AAs = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']

    for seq in seq_types:
        # gets cluster IDs and abundances
        cluster_ids = pd.DataFrame(summary_data[["Clusters", "Percent"]])
        # gets all of the normalized data
        cluster_aa_data = pd.DataFrame(list(summary_data[seq + "_norm"]))
        # adds any missing amino acids with value "0"
        cluster_aa_data[[aa for aa in AAs if aa not in cluster_aa_data]] = 0
        # reorders for easier viewing
        cluster_aa_data = cluster_aa_data[AAs]
        # fills in gaps with 0
        cluster_aa_data = cluster_aa_data.fillna(0)
        # sums all of the AA counts to get the total length of each fragment
        cluster_aa_data["Total Length"] = cluster_aa_data[list(cluster_aa_data.columns)].sum(axis=1)
        # divide each by total length for a percent value
        cluster_aa_data[cluster_aa_data.columns.difference(['Total Length'])] = cluster_aa_data[
            cluster_aa_data.columns.difference(['Total Length'])].div(cluster_aa_data["Total Length"], axis=0)
        # concatenates the IDs and cluster data and reformats for easier reading
        seq_dict[seq] = pd.concat([cluster_ids, cluster_aa_data], axis=1).set_index("Clusters").transpose()
        # save data
        seq_dict[seq].to_excel(writer, sheet_name=seq + "_summary", index=True)

print("Done!")
# endregion
