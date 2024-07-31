# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 | 101 | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 | 8 | +33 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 8 | +33 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 | 101 | +33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem: de-multiplex files and report index-hopping
2. Describe output: 52 files
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

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
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

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
