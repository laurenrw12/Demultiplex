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
    
    a.	Turn in the 4 histograms.
    - link
    - link
    - link
    - link

    b.	What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.

    answer

    c.	How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command(s) you used. CHALLENGE: use a one-line command)

    answer

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

Problem: de-multiplex files and report index-hopping

Input:
- 4 FASTQ files 
```bash
1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz
```
- txt file of 24 matched indexes
```
B1	GTAGCGTA    A5	CGATCGAT    C1	GATCAAGG
B9	AACAGCGA    C9	TAGCCATG    C3	CGGTAATC
B3	CTCTGGAT    C4	TACCGGAT    A11	CTAGCTCA
C7	CACTTCAC    B2	GCTACTCT    A1	ACGATCAG
B7	TATGGCAC    A3	TGTTCCGT    B4	GTCCTAAG
A12	TCGACAAG    C10	TCTTCGAC    A2	ATCATGCG
C2	ATCGTGGT    A10	TCGAGAGT    B8	TCGGATTC
A7	GATCTTGC    B10	AGAGTCCA    A8	AGGATAGC
```
- q score threshold 

Pseudocode:
```python
'''
covert matched index textfile to set

open all 4 input files (using 'a')
  while True:
    read 4 lines (1 record) from all files
    if line == "": break

    reverse compliment index 2
    edit all the headers w/ index 1 seq-index2 seq

    if index 1 & index 2 not in set of 24 known index OR average q score is < q score threshold
      unknown +=1
      write R1 record to unknown_R1.fq
      write R4 record to unknown_R2.fq
      continue
    else:
      if index 1 == index 2:
        matched_dict(index1-index2) += 1
        write R1 record to index_R1.fq
        write R4 record to index_R2.fq
        continue
      else:
        hopped_dict(index1-index2) += 1
        write R1 record to hopped_R1.fq
        write R4 record to hopped_R2.fq
        continue
'''
```

Planned Functions:
```python
def reverse_compliment(sequence: str) -> str:
    '''Takes a DNA sequence and returns the reverse compliment of that sequence.'''
    return rev_complimented_sequence
#Input: CAT
#Expected output: ATG
```
```python
def edit_header(R1_header: str, R4_header: str, R2_sequence: str, R3_sequence: str) -> str, str:
    '''Takes 2 headers and 2 indexes, adds the indexes to both headers, and returns the new headers.'''
    return edited_R1header, edited_R4header
#Input: @K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1, @K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1, NCTTCGAC, NTCGAAGA
#Expected output: @K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1 NCTTCGAC-NTCGAAGA, @K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1 NCTTCGAC-NTCGAAGA
```
```python
def check_qscore(R2_qscores: str, R3_qscores: str, qscore_threshold: int) -> bool:
    '''Takes the qscores of both index files, converts them, averages them, and compares them to the threshold. Returns True or False.'''
    return bool
#Input: EEEEEEE, FFFFFFF, 32
#Expected output: False
```

Output:
52 files
```bash
unknown_R1.fastq
unknown_R2.fastq
hopped_R1.fastq
hopped_R2.fastq
knownindex1_R1.fastq
knownindex1_R2.fastq
knownindex2_R1.fastq
knownindex2_R2.fastq
...
knownindex24_R1.fastq
knownindex24_R2.fastq
```

## The Second 
Looked over classmates puedocode:
- https://github.com/calzamora/Demultiplex/blob/master/Assignment-the-first/pseudo_code_pt2.txt
- https://github.com/lenarayneallen/Demultiplex/blob/master/psuedocode.txt 
- https://github.com/troycho/Demultiplex/blob/master/pseudocode.md

## The Third
Write your code to demultiplex the samples. Be sure to:

- Incorporate feedback from peer code reviews
- Utilize appropriate functions (perhaps you want to `import bioinfo`???)
- Sufficiently comment your code/use docstrings/use type annotations on functions
- Use unit tests on functions/entire algorithm to ensure it works properly
- Create a useful report for the end user of your code
- Use `argparse` to "generalize" your code
- Be mindful of "simple" things you can do to optimize your code
- Follow the specifications laid out in [Assignment the First](../Assignment-the-first#part-2--develop-an-algorithm-to-de-multiplex-the-samples) for the code
    - Unclear? Ask!

Modules that are fair game to import:
- bioinfo
- argparse
- math
- gzip
- numpy
- matplotlib
- itertools

Final work will be submitted on [GitHub in the Assignment-the-Third folder](.). Make sure your folder is well organized and final output is clearly labeled/summarized (a markdown file would be much appreciated!!). Use your code to demultiplex the samples and report:
- Percentage of reads from each sample
- Overall amount of index swapping
- Any figures/any other relevant data your code output
