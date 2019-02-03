#!/usr/bin/env python3

"""

TRS sequence finder for Coronaviruses!

Usage:
    find_trs_mod.py [options] <FASTA> <CSFILE> <OUTPUTFILE>

Options:
    -h, --help                                  Show this little neat help message and exit.
    -r REGION, --region REGION                  Size of extracted region downstream of the TRS. [default: 100]
    -l CS_LENGTH, --length CS_LENGTH            Length of the core sequence. [default: 8]
    -m, --mismatch                              Allow one mismatch in trs sequence. [default: False]
    -t, --withtrs                               Include the trs sequence in extracted sequence. [default: True]

Dependencies:
    BioPython 1.73 or greater
    docopt 0.6.2 or greater

Contact:
    kevin.lamkiewicz{at}uni-jena{dot}de

"""

import csv
import sys

from Bio import SeqIO
from docopt import docopt

def find_all(a_str, sub):
    """


    """
    return([k for k in range(len(a_str)) if a_str[k:k+len(sub)] == sub])    

def hamming(s1, s2):
    """


    """
    if len(s1) == len(s2):
        return int(sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2)))

def find_with_mism(a_str, sub):   #find match with one mismatch allowed
    """


    """
    indexes = []
    k = 0
    while k < len(a_str):
        if hamming(a_str[k:k+len(sub)], sub) == 1 or hamming(a_str[k:k+len(sub)], sub) == 0:  
            indexes.append(k)
        k += 1
    return indexes



args = docopt(__doc__)

file = args['<FASTA>']
csfile = args['<CSFILE>']
outfile = args['<OUTPUTFILE>']
regionSize = int(args['--region'])
cs_length = int(args['--length'])
mismatch = int(args['--mismatch'])
withtrs = int(args['--withtrs'])

sequence = ''
trsSeqs = {}

with open(csfile) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        trsSeqs[row[0]] = row[1]
        #row[0] is the sgmRNA name
        #row[1] is the core sequence
        line_count += 1
    
header = ''
for record in SeqIO.parse(file, "fasta"):
    header = record.id
    sequence = str(record.seq).upper()

if mismatch:
    trsRegion = { geneName : find_all(sequence, trsB) for geneName, trsB in trsSeqs.items() }
else:
    trsRegion = { geneName : find_with_mism(sequence, trsB) for geneName, trsB in trsSeqs.items() }

uniquePositions = set([x for liste in trsRegion.values() for x in liste])

with open(outfile, 'w') as outputStream:
    for trs in uniquePositions:
        if trs - regionSize < 0:
            continue

        if not withtrs:
            subSeq = sequence[trs-regionSize:trs]
            outputStream.write(f">{header}|{trs-regionSize}-{trs}\n{subSeq}\n")
        else:
            subSeq = sequence[trs-regionSize:trs+cs_length] 
            outputStream.write(f">{header}|{trs-regionSize}-{trs+cs_length}\n{subSeq}\n")

