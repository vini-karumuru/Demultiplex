#!/usr/bin/env python

# import libraries
import argparse
import gzip
import numpy as np
import itertools
from itertools import product
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# set global variables to hold inputs
def get_args():
    parser = argparse.ArgumentParser(description="A script that demultiplexes reads from Illumina sequencing files.")
    parser.add_argument("-r1", help="input R1 FASTQ filename", required=True, type=str)
    parser.add_argument("-r2", help="input R2 FASTQ filename", required=True, type=str)
    parser.add_argument("-r3", help="input R3 FASTQ filename", required=True, type=str)
    parser.add_argument("-r4", help="input R4 FASTQ filename", required=True, type=str)
    parser.add_argument("-b", help="text file with barcode list", required=True, type=str)
    parser.add_argument("-o", help="output directory name (must exist already)", required=True, type=str)
    return parser.parse_args()
args = get_args()

# define a function to reverse complement a sequence
def reverse_comp(input_seq: str) -> str:
    # create a dictionary containing complementary base pairs
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    # reverse input sequence
    reverse_seq = input_seq[::-1]
    # initialize an empty string that will contain reverse complemented sequence
    rc_seq = ''
    # loop through bases in reversed input sequence
    for base in reverse_seq:
        # retrieve the base's complement and add this to the reverse complemented sequence string
        rc_seq += comp_dict[base]
    # output the completed reverse complement string
    return rc_seq

# function to write entries to output files
def write_entry(index_seq):
    # retrieve file handles
    fhs = filehandles[index_seq]
    # write R1 entry to <index>_R1.fq
    R1 = entry_dict["R1"]
    fhs[0].write(f"{R1[0]}\n{R1[1]}\n{R1[2]}\n{R1[3]}\n")
    # write R4 entry to <index>_R2.fq
    R4 = entry_dict["R4"]
    fhs[1].write(f"{R4[0]}\n{R4[1]}\n{R4[2]}\n{R4[3]}\n")
    # returns nothing

# parse through barcodes file to get list of barcodes
indexes = []
with open(args.b, 'r') as fh:
    for ind,line in enumerate(fh):
        if ind > 0:
            index_seq = line.strip('\n').split('\t')[4]
            indexes.append(index_seq)

# initialize dictionary where keys are combination of valid indexes + values are counts of index combination
ind_combo_dict = {}
# get all possible permutations (with replacement) of indexes
index_permutations = product(indexes, repeat=2)
for perm in index_permutations:
    # combine indexes with "-" in between
    index_combo = f"{perm[0]}-{perm[1]}"
    # initialize count of index combination to 0
    ind_combo_dict[index_combo] = 0

# initialize summary counters
correct_dual = 0
index_hopped = 0
unknown_indexes = 0

# get number of entries/paired reads in sequencing dataset
with gzip.open(args.r2, "r") as fh:
    # loop through lines in file
    for lin_num, line in enumerate(fh):
        pass
num_entries = (lin_num + 1)/4

file_prefixes = indexes + ['undetermined', 'index_hopped']
# create a dictionary to store file handles: 
## key is name of index
## value is a list of corresponding R1 + R2 file handles
filehandles = {}
for prefix in file_prefixes:
    # open R1 + R2 output FASTQ files for writing
    file_handle1 = open(f"{args.o}/{prefix}_R1.fq", "w")
    file_handle2 = open(f"{args.o}/{prefix}_R2.fq", "w")
    # save filehandles to dictionary
    filehandles[prefix] = [file_handle1, file_handle2]

# open FASTQ files for reading
with gzip.open(args.r1, "r") as r1, gzip.open(args.r2, "r") as r2, gzip.open(args.r3, "r") as r3, gzip.open(args.r4, "r") as r4:
    entries = ["R1", "R2", "R3", "R4"]

    # loop through each entry
    for entry in np.arange(num_entries):
        # initialize/empty dictionary that will contain current entry
        entry_dict = {"R1": [], "R2": [], "R3": [], "R4": []}
        # loop through 4 input file handles
        for ind, fh in enumerate([r1, r2, r3, r4]):
            # slice out 4 lines (1 entry) & looping through these lines
            for line in itertools.islice(fh, 4):
                # convert line (currently a bytes object) to a string & strip new line character
                line = line.decode('utf-8').strip('\n')
                # append line to corresponding entry
                entry_dict[entries[ind]].append(line)

        # replace sequence of R3 with its reverse complement
        R3_rc = reverse_comp(entry_dict["R3"][1])
        entry_dict["R3"][1] = R3_rc

        R2_ind = entry_dict["R2"][1]
        # combine R2 and R3 indexes into 1 string separated by "-"
        combined_ind = f"{R2_ind}-{R3_rc}"
        # add R2 and R3 indexes to header lines of R1 and R4
        entry_dict["R1"][0] += f" {combined_ind}"
        entry_dict["R4"][0] += f" {combined_ind}"

        # sort any reads containing "N" into undetermined FASTQ files
        if "N" in combined_ind:
            write_entry("undetermined")
            unknown_indexes += 1
            # skip to next entry
            continue

        # if both indexes are not valid indexes
        if set([R2_ind, R3_rc]).issubset(set(indexes)) == False:
            # sort entry into undetermined FASTQ files
            write_entry("undetermined")
            unknown_indexes += 1
            # skip to next entry
            continue

        ind_combo_dict[combined_ind] += 1
        # if both indexes are the same valid index
        if R2_ind == R3_rc:
            # sort entry into FASTQ files for that index
            write_entry(R2_ind)
            correct_dual += 1
        # if both indexes are valid but different
        else:
            # sort entry into index hopped FASTQ files
            write_entry("index_hopped")
            index_hopped += 1

# closing all files
for file_prefix in filehandles:
    fh_pair = filehandles[file_prefix]
    fh_pair[0].close()
    fh_pair[1].close()

# write summary counts & percentages to .txt output file
with open(f"{args.o}/read_count_summary.txt", 'w') as sfh:
    sfh.write(f"Number of read-pairs with properly matched indexes (successfully demultiplexed): {correct_dual} ({round(correct_dual/num_entries*100, 2)}%)\n")
    sfh.write(f"Number of read pairs with index-hopping observed: {index_hopped} ({round(index_hopped/num_entries*100, 2)}%)\n")
    sfh.write(f"Number of read-pairs with unknown index(es): {unknown_indexes} ({round(unknown_indexes/num_entries*100, 2)}%)\n")

# create matrix of counts
count_array = np.zeros((len(indexes),len(indexes)), dtype = int)
for i, index1 in enumerate(indexes):
    for j, index2 in enumerate(indexes):
        count = ind_combo_dict[f"{index1}-{index2}"]
        count_array[i, j] = count

# writing matrix to markdown file
with open(f"{args.o}/index_combination_counts.md", 'w') as md:
    markdown_text = ""
    header_row = "| | " + " | ".join(indexes) + " |\n"
    separator_line =  "| " + " | ".join(["---"] * (len(indexes)+1)) + " |\n"
    markdown_text += header_row + separator_line
    for i, index in enumerate(indexes):
        index_counts = list(count_array[i,:].astype(str))
        index_row = "| " + index + " | " + " | ".join(index_counts) + " |\n"
        markdown_text += index_row
    md.write(markdown_text)

# creating a heatmap of matrix
fig, ax = plt.subplots()
im = ax.imshow(count_array, cmap='plasma', norm=LogNorm())
# adding a color bar
cbar = plt.colorbar(im)
cbar.set_label('Log (Count)')
# labeling title & axes
ax.set_title('Heatmap of Log Index Combination Counts')
ax.set_ylabel('Index 1')
ax.set_xlabel('Index 2')
ax.set_xticks(range(len(indexes)), labels=indexes, rotation=45, ha='right', rotation_mode="anchor")
ax.set_yticks(range(len(indexes)), labels=indexes)
# adding number labels
# for i in range(len(indexes)):
#     for j in range(len(indexes)):
#         text = ax.text(j, i, count_array[i, j], ha='center', va='center', color='w')
fig.tight_layout()
plt.savefig(f"{args.o}/index_combination_heatmap.png") 
plt.close()

# creating a bar chart of reads per sample
index_counts = [ind_combo_dict[f"{index}-{index}"] for index in indexes]
plt.bar(indexes, index_counts)
plt.xlabel("Index")
plt.ylabel("Number of Paired Reads")
plt.title("Number of Paired Reads per Sample/Index")
plt.xticks(rotation=45)
fig.tight_layout()
plt.savefig(f"{args.o}/reads_per_sample_barchart.png")
plt.close()
