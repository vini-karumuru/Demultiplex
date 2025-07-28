#!/usr/bin/env python

# import libraries
import argparse
import gzip
import bioinfo
import matplotlib.pyplot as plt

# set global variables to hold inputs
def get_args():
    parser = argparse.ArgumentParser(description="A script that takes in a FASTQ file and produces a plot of the per base quality score distribution.")
    parser.add_argument("-i", help="input FASTA filename", required=True, type=str)
    parser.add_argument("-o", help="output plot filename", required=True, type=str)
    return parser.parse_args()
args = get_args()

# opening & reading gzipped input FASTQ file
with gzip.open(args.i, 'rb') as fh:

    # determine read length
    for line_ind, line in enumerate(fh):
        if line_ind == 3:
            # converting line (currently a bytes object) to a string and stripping new line character
            line = line.decode('utf-8').strip('\n')
            read_len = len(line)
            # stop looping through file once read length is calculated
            break
        
    # creating a list to hold sum of quality scores at each position, initializing each item to 0.0
    sum_list = [0.0 for i in range(read_len)]

    # sum up quality scores at each position
    for line_ind, line in enumerate(fh):
        # only working with quality score lines
        if (line_ind + 1) % 4 == 0:
            # converting line (currently a bytes object) to a string and stripping new line character
            line = line.decode('utf-8').strip('\n')
            # looping through quality scores in line
            for position_ind, letter in enumerate(line):
                # decoding the quality score
                qs = bioinfo.convert_phred(letter)
                # adding quality score to the running sum of quality scores at that position in read
                sum_list[position_ind] += qs

    # determining number of records in FASTQ file
    num_records = int(line_ind + 1)/4

    # averaging quality scores at each position by dividing the summed quality score by the number of records
    avg_list = [summed_qs/num_records for summed_qs in sum_list]
    
# create a plot of quality score by position
x_values = range(read_len)
plt.figure(figsize=(14, 6))
plt.plot(x_values, avg_list)
plt.xlabel('Position')
plt.ylabel('Average Quality Score')
plt.title('Average Quality Score by Position in Read')
plt.savefig(args.o)