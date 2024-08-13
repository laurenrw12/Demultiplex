# Demultiplexing Lab Notebook

## The First

### Part 1 – Quality Score Distribution per-nucleotide
1.	Perform some initial data exploration! Record any bash commands you used inside a lab notebook (submit to this repo!).

Look at each of the file using zcat:
```
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1
GNCTGGCATTCCCAGAGACATCAGTACCCAGTTGGTTCAGACAGTTCCTCTATTGGTTGACAAGGTCTTCATTTCTAGTGATATCAACACGGTGTCTACAA
+
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 2:N:0:1
NCTTCGAC
+
#AA<FJJJ

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 3:N:0:1
NTCGAAGA
+
#AAAAJJF

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1
NTTTTGATTTACCTTTCAGCCAATGAGAAGGCCGTTCATGCAGACTTTTTTAATGATTTTGAAGACCTTTTTGATGATGATGATGTCCAGTGAGGCCTCCC
+
#AAFAFJJ-----F---7-<FA-F<AFFA-JJJ77<FJFJFJJJJJJJJJJAFJFFAJJJJJJJJFJF7-AFFJJ7F7JFJJFJ7FFF--A<A7<-A-7--
```





1.	Initial Data Exploration Continued:

    a.	Determine which files contain the indexes, and which contain the paired end reads containing the biological data of interest. Create a table and label each file with either read1, read2, index1, or index2.
	
    b.	Determine the length of the reads in each file.
	
    c.	Determine the phred encoding for these data.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 | 101 | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 | 8 | +33 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 8 | +33 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 | 101 | +33 |

Phred Encoding: <img width="831" alt="Screenshot 2024-07-30 at 10 39 05 AM" src="https://github.com/user-attachments/assets/074fca78-a26f-4807-ac8e-11df826817ec">

2.	Generate a per base distribution of quality scores for read1, read2, index1, and index2. Average the quality scores at each position for all reads and generate a per nucleotide mean distribution **as you did in part 1 of PS4 in Bi621**. (NOTE! Do NOT use the 2D array strategy from PS9 - you WILL run out of memory!)

    [ANSWERS.md](Answers.md)

### Part 2 – Develop an algorithm to de-multiplex the samples
Write up a strategy (**NOT A SCRIPT**) for writing an algorithm to de-multiplex files and reporting index-hopping. That is, given four input FASTQ files (2 with biological reads, 2 with index reads) and the 24 known indexes above, demultiplex reads by index-pair, outputting:

- one R1 FASTQ file and one R2 FASTQ file **per** matching index-pair, 
- another two FASTQ files for non-matching index-pairs (index-hopping), and 
- two additional FASTQ files when one or both index reads are unknown or low quality (do not match the 24 known indexes [this includes indexes with 'N's in them] or do not meet a quality score cutoff)
    
Add the sequence of the index-pair to the header of BOTH reads in all of your FASTQ files for all categories (e.g. add “AAAAAAAA-CCCCCCCC” to the end of headers of every read pair that had an index1 of AAAAAAAA and an index2 of CCCCCCCC; this pair of reads would be in the unknown category as one or both of these indexes do not match the 24 known indexes).

Additionally, your algorithm should report: 
- the number of read-pairs with properly matched indexes (per index-pair), 
- the number of read pairs with index-hopping observed, and
- the number of read-pairs with unknown index(es).

You should strive to report values for each possible pair of indexes (both swapped and dual matched). **You should not write any code for this portion of the assignment**. 

[pseudocode.md](pseudocode.md)

## The Second 
Looked over classmates puedocode:
- https://github.com/calzamora/Demultiplex/blob/master/Assignment-the-first/pseudo_code_pt2.txt
- https://github.com/lenarayneallen/Demultiplex/blob/master/psuedocode.txt 
- https://github.com/troycho/Demultiplex/blob/master/pseudocode.md

## The Third
Write your code to demultiplex the samples. Be sure to:

- Incorporate feedback from peer code reviews
    - one big if statement with elif statements instead of many nested if statements

- Utilize appropriate functions 
    - import bioinfo, argparse, gzip 

- Sufficiently comment your code/use docstrings/use type annotations on functions
- Use unit tests on functions/entire algorithm to ensure it works properly
- Create a useful report for the end user of your code
- Use `argparse` to "generalize" your code
- Be mindful of "simple" things you can do to optimize your code
- Follow the specifications laid out in [Assignment the First](../Assignment-the-first#part-2--develop-an-algorithm-to-de-multiplex-the-samples) for the code. Unclear? Ask!

Final work will be submitted on [GitHub in the Assignment-the-Third folder](.). Make sure your folder is well organized and final output is clearly labeled/summarized (a markdown file would be much appreciated!!). Use your code to demultiplex the samples and report:
- Percentage of reads from each sample
- Overall amount of index swapping
- Any figures/any other relevant data your code output

### 08/06/24 Notes
- 

### Final Output Statistics
#### Output for Q-score Threshold of 20:
Total Number of Reads: 363246735

Amount of Unknown Reads: 31828861\
Percentage of Unknown Reads: 8.762325420488638%

Amount of Matched Reads: 330738415\
Percentage of Matched Reads: 91.05062293264659%

Amount of Hopped Reads: 679459\
Percentage of Hopped Reads: 0.1870516468647681%

#### Output for Q-score Threshold of 30:
Total Number of Reads: 363246735

Amount of Unknown Reads: 57748853\
Percentage of Unknown Reads: 15.897968910856141%

Amount of Matched Reads: 304980270\
Percentage of Matched Reads: 83.95953510772782%

Amount of Hopped Reads: 517612\
Percentage of Hopped Reads: 0.1424959814160477%