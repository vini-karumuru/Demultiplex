#!/usr/bin/env python

import bioinfo

indexes = []
with open('/projects/bgmp/shared/2017_sequencing/indexes.txt', 'r') as fh:
    for ind,line in enumerate(fh):
        if ind > 0:
            index_seq = line.strip('\n').split('\t')[4]
            indexes.append(index_seq)

def calc_hamming(seq1, seq2):
    hamming_dist = 0
    for index in range(len(seq1)):
        if seq1[index] != seq2[index]:
            hamming_dist += 1
    return hamming_dist

hamming_distances = []
for ind1, index1 in enumerate(indexes):
    for ind2, index2 in enumerate(indexes):
        if ind1 != ind2:
            hamming_distance = calc_hamming(index1, index2)
            hamming_distances.append(hamming_distance)

print(f"The lowest Hamming distance is {min(hamming_distances)}")
print(f"The median Hamming distance is {bioinfo.calc_median(sorted(hamming_distances))}")
print(f"The highest Hamming distance is {max(hamming_distances)}")
