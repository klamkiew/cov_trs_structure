#!/usr/bin/env python3

"""

TRS sequence finder for Coronaviruses!

Usage:
    find_trs_mod.py [options] <FASTA> <CSFILE> <OUTPUTFILE>

Options:
    -h, --help                                  Show this little neat help message and exit.
    -r REGION, --region REGION                  Size of extracted region downstream of the TRS. [default: 150]
    -m, --mismatch                              Allow one mismatch in trs sequence. [default: False]
    -t, --withtrs                               Include the trs sequence in extracted sequence. [default: False]

Dependencies:
    BioPython 1.73 or greater
    docopt 0.6.2 or greater

Contact:
    kevin.lamkiewicz{at}uni-jena{dot}de

"""

###################################

import csv
import sys

from Bio import SeqIO
from docopt import docopt

###################################

def find_all(haystack, needle):
    """
    Small helper-function to find all occurrences of a substring
    within a string. Returns the starting index in a similar fashion 
    as the built-in str.find() function.

    Parameters:
    haystack -- String that is scanned for the needle
    needle -- the pattern / substring of interest 

    Return:
    List of all starting indices of needle in haystack.

    """
    return([k for k in range(len(haystack)) if haystack[k:k+len(needle)] == needle])    

def hamming(s1, s2):
    """
    Returns the Hamming distance of two input strings.
    Careful, only strings of equal length are allowed.

    Parameters:
    s1 -- First string for distance calculation
    s2 -- Second string for distance calculation

    Return:
    Integer value containing the hamming distance.

    """
    if len(s1) == len(s2):
        return int(sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2)))

def find_with_mism(haystack, needle):   #find match with one mismatch allowed
    """
    Finds all occurences of needle in haystack even if there is one mismatch
    between the pattern and the text. Uses the hamming() function.

    Parameters:
    haystack -- String that is scanned for the needle
    needle -- the pattern / substring of interest

    Return:
    List of all starting indices of needle in haystack with at most one mismatch.

    """
    indexes = []
    k = 0
    while k < len(haystack)-k:
        if hamming(haystack[k:k+len(needle)], needle) <= 1:  
            indexes.append(k)
        k += 1
    return indexes

###################################

if __name__ == '__main__':
    """
    Main method.
    """
    ###################################
    # Argument parsing
    args = docopt(__doc__)

    file = args['<FASTA>']
    csfile = args['<CSFILE>']
    outfile = args['<OUTPUTFILE>']
    regionSize = int(args['--region'])
    #cs_length = int(args['--length'])
    mismatch = int(args['--mismatch'])
    withtrs = args['--withtrs']

    ###################################
    # variable declaration

    sequence = ''
    trsSeqs = {}

    ###################################
    # get CS information from csv

    with open(csfile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if row:
            # row[0] is the sg mRNA name
            # row[1] is the core sequence
            # row[2] is the starting position within the genome
                trsSeqs[row[0]] = (row[1].upper().replace('U','T'), row[2])
            
            
    
    ###################################
    # read in reference genome

    header = ''
    for record in SeqIO.parse(file, "fasta"):
        header = record.id
        sequence = str(record.seq).upper()

    ###################################
    # get all genomic positions
    # for each TRS sequences from the
    # reference genome 

    if mismatch:
        trsRegion = { geneName : find_with_mism(sequence, trsB[0]) for geneName, trsB in trsSeqs.items() }
    else:
        trsRegion = { geneName : find_all(sequence, trsB[0]) for geneName, trsB in trsSeqs.items() }

    ###################################
    # write down results

    with open(outfile, 'w') as outputStream:
        for geneName, trsB in trsSeqs.items():    
            try:
                # just consider the positions
                # that are covered by the literature to be TRS-B
                pos = int(trsB[1])
                cs_length = len(trsB[0])
                canonicalTRS = [x for x in trsRegion[geneName] if x in range(pos-40, pos+40)]

                for trs in canonicalTRS:
                    if trs - regionSize < 0:
                        continue
                    if not withtrs:
                        subSeq = sequence[trs-regionSize:trs]
                        coordinates = f"{trs-regionSize}-{trs}"
                    else:
                        subSeq = sequence[trs-regionSize:trs+cs_length] 
                        coordinates = f"{trs-regionSize}-{trs+cs_length}"
                    outputStream.write(f">{header}|{coordinates}|Gene_{geneName}\n{subSeq}\n")
            except ValueError:
                continue
    exit(0)
