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

  

