#!/usr/bin/env python
import argparse
import bioinfo
from matplotlib import pyplot as plt 
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", help="Name of input file", required=True, type=str)
    parser.add_argument("-l", help="length of sequences", required=True, type=int)
    parser.add_argument("-o", help="Name of output file", required=True, type=str)
    return parser.parse_args() 

args = get_args()
f = args.f
l = args.l
o = args.o

#2.	Generate a per base distribution of quality scores for read1, read2, index1, and index2. 
# Average the quality scores at each position for all reads and generate a per nucleotide mean distribution 
# **as you did in part 1 of PS4 in Bi621**. 
# (NOTE! Do NOT use the 2D array strategy from PS9 - you WILL run out of memory!)

mean_qscores: list = []
for i in range (0,l):
        mean_qscores.append(0)

with gzip.open(f,"rt") as f:
    i = 0
    for line in f:
        i+=1
        #line = str(line)
        line = line.strip('\n')
        if (i-1)%4 == 3:
            for j in range(len(line)): 
                mean_qscores[j] += bioinfo.convert_phred(line[j])
    lines = i

num_records = int(lines/4)

for i, sum in enumerate(mean_qscores):
    average = sum/(num_records)
    mean_qscores[i] = average

#THIS IS MY PRINT CODE
#print("# Base Pa#ir\tMean Quality Score")
#for x in range(len(mean_qscores)):
#    print(x, mean_qscores[x], sep="\t")

#THIS IS MY GRAPH
x = range(0,l)
y = mean_qscores

plt.bar(x,y) 

plt.xlabel("Base Pair Number")
plt.ylabel("Mean Quality Score")
plt.title("Mean Quality Score by Base Pair")
plt.savefig(o)