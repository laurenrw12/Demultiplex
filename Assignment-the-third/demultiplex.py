#!/usr/bin/env python
import argparse
import bioinfo
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f1", help="Name of 1st input file", required=True, type=str)
    parser.add_argument("-f2", help="Name of 2nd input file", required=True, type=str)
    parser.add_argument("-f3", help="Name of 3rd input file", required=True, type=str)
    parser.add_argument("-f4", help="Name of 4th input file", required=True, type=str)
    parser.add_argument("-q", help="Q-Score Threshold", required=True, type=int)
    parser.add_argument("-m", help="Name of matched index txt file", required=True, type=str)
    return parser.parse_args() 

args = get_args()
f1 = args.f1
f2 = args.f2
f3 = args.f3
f4 = args.f4
q = args.q
m = args.m

def reverse_complement(seq: str) -> str:
    '''Takes a DNA sequence and returns the reverse complement of that sequence.'''
    #initialize 
    pairs = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
    rev_comp_seq = ""

    #reverse the string
    rev_seq = seq[::-1]

    #complement the string & add to rev_comp_seq 
    for base in rev_seq:
        rev_comp_seq += pairs[base]
    return rev_comp_seq 
#print(reverse_complement("CAT"))
#Input: CAT
#Expected output: ATG

def check_quality(R2_qscores: str, R3_qscores: str, qscore_threshold: int) -> bool:
    '''Takes the qscores of both index files, converts them, averages them, and compares them to the threshold. Returns True (quality) or False (not quality).'''
    convertedR2 = bioinfo.qual_score(R2_qscores)
    convertedR3 = bioinfo.qual_score(R3_qscores)
    if convertedR2 < qscore_threshold or convertedR3 < qscore_threshold:
        return False
    else:
        return True
#print(check_quality("EEEEEEE", "FFFFFFF", 32))
#Input: EEEEEEE, FFFFFFF, 32
#Expected output: False

def write_to_files(key:str) -> None:
    '''takes a key from write_dict'''
    #write R1 record to key_R1.fq
    outfileR1 = write_dict[key][0]
    outfileR1.write(f'{f1header}\n{f1seq}\n{f1plus}\n{f1qscore}\n')
    #write R4 record to key_R2.fq
    outfileR2 = write_dict[key][1]
    outfileR2.write(f'{f4header}\n{f4seq}\n{f4plus}\n{f4qscore}\n')
#Input: "hopped"
#Expected output:

#initial variables
matched_indexes = set()
matched_dict:dict[str,int] = {}
hopped_dict:dict[str,int] = {}
write_dict:dict = {}
unknown = 0

#covert matched index textfile to set
with open(m, "r") as matched:
    for line in matched:
        line = line.strip("\n")
        line = line.split("\t")
        matched_indexes.add(line[1])
        matched_indexes.add(line[3])
        matched_indexes.add(line[5])

#open all 4 input files 
with gzip.open(f1,"rt") as f1, gzip.open(f2,"rt") as f2, gzip.open(f3,"rt") as f3, gzip.open(f4,"rt") as f4: 

    #open all 52 files for writing
    for index in matched_indexes:
        fh1 = open(f'./output_{q}/{index}_R1.fastq', "a")
        fh2 = open(f'./output_{q}/{index}_R2.fastq', "a")
        write_dict[index] = (fh1, fh2)

    fh1 = open(f'./output_{q}/hoppped_R1.fastq', "a")
    fh2 = open(f'./output_{q}/hoppped_R2.fastq', "a")
    write_dict["hopped"] = (fh1, fh2)

    fh3 = open(f'./output_{q}/unknown_R1.fastq', "a")
    fh4 = open(f'./output_{q}/unknown_R2.fastq', "a")
    write_dict["unknown"] = (fh3, fh4)

    while True:
        #read 4 lines (1 record) from all files
        f1header, f1seq, f1plus, f1qscore = f1.readline().strip(), f1.readline().strip(), f1.readline().strip(), f1.readline().strip()
        f2header, f2seq, f2plus, f2qscore = f2.readline().strip(), f2.readline().strip(), f2.readline().strip(), f2.readline().strip()
        f3header, f3seq, f3plus, f3qscore = f3.readline().strip(), f3.readline().strip(), f3.readline().strip(), f3.readline().strip()
        f4header, f4seq, f4plus, f4qscore = f4.readline().strip(), f4.readline().strip(), f4.readline().strip(), f4.readline().strip()

        #end while True loop
        if f1header == "": 
            break

        #reverse compliment index 2
        f3seq = reverse_complement(f3seq) # type: ignore

        #create combo index of f2seq_f3seq
        f2seq_f3seq = f'{f2seq}-{f3seq}'

        #edit  the headers w/ index 1 seq-index2 seq
        f1header, f4header = f1header + " " + f2seq_f3seq, f4header + " " + f2seq_f3seq # type: ignore

        #if (T & T) or not F = if T or T = T         sequences not in matched & not quality -> UNKNOWN
        #if (T & T) or not T = if T or F = T         sequences not in matched & quality -> UNKNOWN
        #if (F & F) or not F = if F or T = T         sequences in matched, & not quality -> UNKNOWN
        #if (F & F) or not T = if F or F = F         sequences in matched & quality -> GOOD, KEEP GOING
        if (f2seq not in matched_indexes or f3seq not in matched_indexes) or not check_quality(f2qscore, f3qscore, q): # type: ignore
            unknown +=1
            write_to_files("unknown")
            continue

        elif f2seq == f3seq:
            if f2seq_f3seq in matched_dict:
                matched_dict[f2seq_f3seq] += 1
            else: 
                matched_dict[f2seq_f3seq] = 1 

            write_to_files(f2seq) # type: ignore
            continue

        else:
            if f2seq_f3seq in hopped_dict:
                hopped_dict[f2seq_f3seq] += 1
            else: 
                hopped_dict[f2seq_f3seq] = 1

            write_to_files("hopped")

    #close all the open files
    for filename in write_dict:
        write_dict[filename][0].close()
        write_dict[filename][1].close()

#calculate total # of matched reads
matched = 0
for key in matched_dict:
    matched += matched_dict[key]

#calculate total # of matched reads
hopped = 0
for key in hopped_dict:
    hopped += hopped_dict[key]

#total # of reads
totalreads = matched + hopped + unknown 

#report statistics
print(f"Total Number of Reads: {totalreads}\n")
print(f"Amount of Unknown Reads: {unknown}\tPercentage of Unknown Reads: {(unknown/totalreads)*100}%\n")
print(f"Amount of Matched Reads: {matched}\tPercentage of Matched Reads: {(matched/totalreads)*100}%\n")
print(f"Amount of Hopped Reads: {hopped}\tPercentage of Hopped Reads: {(hopped/totalreads)*100}%\n")

print("Amount of Reads per Sample")
for key in matched_dict:
    print(f"Sample: {key}\tAmount of Reads: {matched_dict[key]}\tPercentage of Reads: {(matched_dict[key]/totalreads)*100}%")