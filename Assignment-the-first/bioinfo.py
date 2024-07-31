#!/usr/bin/env python

# Author: Lauren Williams
# Email: lwil@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.2"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set('ACGNTacgnt')
RNA_bases = set('ACGNUacgnu')

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    '''Calculates the average quality score of a whole phred string'''
    phred_score_list = list(phred_score)
    sum = 0 
    for character in range(0,len(phred_score_list)):
        sum += convert_phred(phred_score_list[character]) 
    return sum/len(phred_score_list)

def validate_base_seq(seq:str, RNAflag:bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(seq:str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(seq) #String contains invalid characters - are you sure you used a DNA sequence?
    seq = seq.upper()
    return (seq.count("G")+seq.count("C"))/len(seq)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    length = len(lst)
    if length % 2 == 0:
        index = length //2 
        return ((lst[index-1] + lst[index]) /2)
    else: 
        index = length // 2
        return (lst[index])
        

def oneline_fasta(inputfile:str) -> str:
    '''Given a file, return the name of a new file with all the sequence strings on ONE line. i.e. two lines per record '''
    with open(f'oneline_{inputfile}',"w") as outf:
        with open(inputfile, "r") as inf:
            linenum = 0
            while True:
                line = inf.readline().strip()
                if line == "":
                    break

                if line.startswith(">") and linenum == 0:
                    outf.write(f'{line}\n')
                    linenum += 1
                elif line.startswith(">") and linenum != 0:
                    outf.write(f'\n{line}\n')
                    linenum += 1
                else: 
                    outf.write(line)
                    linenum += 1
    return f'oneline_{inputfile}'

if __name__ == "__main__":
    # write tests for functions above
    # tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")