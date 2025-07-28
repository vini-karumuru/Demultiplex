#!/usr/bin/env python

# Author: Vini Karumuru vini@uoregon.edu

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.'''

__version__ = "0.1"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning


def convert_phred(letter: str) -> int:
    '''
    This function takes in a single character (must be a string) that is an ASCII-encoded phred+33 quality score. 
    It converts the character into a phred score and returns this score as an integer.
    '''
    assert len(letter) == 1, "Input must be a single character."
    return ord(letter) - 33

def qual_score(phred_scores: str) -> float:
    '''
    This function takes in a string of ASCII-encoded phred+33 quality scores. 
    It will convert the scores to Illumina quality scores and return the average quality score as a float.
    '''
    summed_score = 0
    for char in phred_scores:
        summed_score += convert_phred(char)
    return summed_score/len(phred_scores)

DNA_bases = set("ATCGNatcgn")
RNA_bases = set("AUCGNaucgn")
def validate_base_seq(seq: str, RNAflag: bool = False):
    '''
    This function takes in
    (1) a nucleotide sequence (string) and 
    (2) a boolean value (True or False) to signify if the sequence is RNA. Default is False, which would signify DNA input.
    Returns True if the string is composed of only As, Ts (or Us if RNAflag is True), Gs, and Cs. 
    False is returned otherwise.
    Case insensitive.
    '''
    seq = seq.upper()
    return set(seq) <= (RNA_bases if RNAflag else DNA_bases)

    

def gc_content(seq: str) -> bool:
    '''
    This function takes in a DNA or RNA sequence (string) and returns the GC content of the sequence
    (the proportion of bases that are G or C) as a decimal (between 0 and 1). 
    Case insensitive.
    '''
    assert validate_base_seq(seq, True) or validate_base_seq(seq, False), "Not a valid nucleotide sequence. Sequence must only be comprised of As, Cs, Ts (or Us if RNA), and/or Gs."
    seq = seq.upper()
    GC_count = 0
    for char in seq:
        if (char == 'G') or (char == 'C'):
            GC_count += 1
    return GC_count/len(seq)

def calc_median(input_list: list) -> float:
    '''Given a sorted list, this function returns the median value of the list as a float.'''
    assert sorted(input_list) == input_list, "Input list is not sorted."
    if len(input_list) % 2 == 1:
        median_position = len(input_list) // 2
        return input_list[median_position]
    else:
        median_position_1 = len(input_list) // 2
        median_position_2 = median_position_1 -1
        median_value = (input_list[median_position_1] + input_list[median_position_2])/2
        return median_value

def oneline_fasta(input_fasta: str, output_filename: str):
    '''
    This function takes in a FASTA file and an output filename.
    The function combines the sequence lines for any entry in the FASTA file into a single line and 
    outputs a FASTQ file in which each entry is exactly 2 lines (the header and the sequence).
    '''
    assert ('.fq' in input_fasta) or ('.fasta' in input_fasta), "Input file must be a fasta (.fa or .fasta) file."
    assert ('.fq' in output_filename) or ('.fasta' in output_filename), "Output filename must contain '.fa' or '.fasta.'"
    with open(output_filename, 'w') as o_fh:
        with open(input_fasta, 'r') as i_fh:
            sequence_lines: str = ""
            for line in i_fh:
                if ">" in line:
                    if sequence_lines == "":
                        o_fh.write(line)
                    else:
                        sequence_lines += "\n"
                        o_fh.write(sequence_lines)
                        sequence_lines = ""
                        o_fh.write(line)
                else:
                    line = line.strip("\n")
                    sequence_lines += line
            o_fh.write(sequence_lines)

if __name__ == "__main__":
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job!")

    assert qual_score("I") == 40, "wrong average phred score for 'I'"
    assert qual_score("I$") == 43/2, "wrong average phred score for 'I$'"
    assert qual_score("IAD") == 107/3, "wrong average phred score for 'IAD'"
    assert qual_score("ACB&6") == 25, "wrong average phred score for 'ACB&6'"
    assert qual_score("BCCCCCCB") == 33.75, "wrong average phred score for 'BCCCCCCB'"
    print("Your qual_score function is working! Nice job!")

    assert validate_base_seq("T", False) == True, "your function does not work with single letter sequences"
    assert validate_base_seq("UaseTfwr2A") == False, "'UaseTfwr2A' is not a valid DNA sequence but your function said it was"
    assert validate_base_seq("ATCTT", False) == True, "'ATCTT' is a valid DNA sequence but your function said it was not"
    assert validate_base_seq("AUUCG", True) == True, "'AUUCG' is a valid RNA sequence but your function said it was not"
    assert validate_base_seq("AUU C", True) == False, "'AUU C' is not a valid RNA sequence because of the space but your function said it was"
    assert validate_base_seq("aaaaat", False) == True, "lowercase sequences are valid base sequences"
    assert validate_base_seq("aGaCa") == True, "lowercase sequences are valid base sequences"
    print("Your validate_base_seq function is working! Nice job!")

    assert gc_content("G") == 1, "your function does not work with single letter sequences"
    assert gc_content("GU") == 0.5, "incorrect GC content for 'GU'"
    assert gc_content("TATTTTA") == 0, "incorrect GC content for 'TATTTTA'"
    assert gc_content("GCGCGCGCGGGGGGCCCCC") == 1, "incorrect GC content for 'GCGCGCGCGGGGGGCCCCC'"
    assert gc_content("TATTGCGGAGCTC") == 7/13, "incorrect GC content for 'TATTGCGGAGCTC'"
    print("Your gc_content function is working! Nice job!")

    assert calc_median([1]) == 1, "incorrect median of [1]"
    assert calc_median([23,46,990]) == 46, "incorrect median of [23,46,990]"
    assert calc_median([350,400,453,873]) == 853/2, "incorrect median of [350,400,453,873]"
    assert calc_median([2,2,2.5,2.5,2.78,3,3.5]) == 2.5, "incorrect median of [2,2,2.5,2.5,2.78,3,3.5]"
    assert calc_median([134.1423, 413.5342, 634.8452, 5234.34]) == 524.1897, "incorrect median of [134.1423, 413.5342, 634.8452, 5234.34]"
    print("Your calc_median function is working! Nice job!")
