#!/usr/bin/env python
import argparse
import bioinfo
from matplotlib import pyplot as plt 
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-file", help="Name of input file", required=True, type=str)
    return parser.parse_args() 

args = get_args()
file = args.file

#2.	Generate a per base distribution of quality scores for read1, read2, index1, and index2. 
# Average the quality scores at each position for all reads and generate a per nucleotide mean distribution 
# **as you did in part 1 of PS4 in Bi621**. 
# (NOTE! Do NOT use the 2D array strategy from PS9 - you WILL run out of memory!)

mean_qscores: list = []
for i in range (0,8):
        mean_qscores.append(0.0)

with gzip.open(file,"rt") as f:
    i = 0
    for line in f:
        i+=1
        line = str(line)
        line = line.strip('\n')
        if (i-1)%4 == 3:
            for j, character in enumerate(line): 
                mean_qscores[j] += bioinfo.convert_phred(character)
    lines = i

num_records = int(lines/4)

for i, sum in enumerate(mean_qscores):
    average = sum/(num_records)
    mean_qscores[i] = average

#THIS IS MY PRINT CODE
#print("# Base Pa#ir\tMean Quality Score")
#for x in range(len(mean_qscores)):
#    print(x,"\t",mean_qscores[x])

#THIS IS MY GRAPH
x = range(0,101)
y = mean_qscores

plt.bar(x,y) 
plt.xlabel("Base Pair Number")
plt.ylabel("Mean Quality Score")
plt.title("Mean Quality Score by Base Pair")