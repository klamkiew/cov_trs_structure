#!/usr/bin/env python3

"""

TRS-L matching to genome to find putative hybridization subsequences.

Usage:
    leader_hybridization.py [options] <GENOME_FASTA> <OUTPUTFILE>

Options:
    -h, --help                          Show this little neat help message and exit.
    -r FLANKING, --region FLANKING      Size of extracted region up and downstream of the TRS_L. [default: 5]
    -l CS_LENGTH, --length CS_LENGTH    Length of the core sequence. [default: 8]
"""

from Bio import SeqIO
import RNA
import sys
import csv
import subprocess
from docopt import docopt

def find_all(a_str, sub):
    return([k for k in range(len(a_str)) if a_str[k:k+len(sub)] == sub])    

file = args['<GENOME_FASTA>']
outfile = args['<OUTPUTFILE>']
regionSize = int(args['--region'])
cs_length = int(args['--length'])

trsL = 'CGTTTAGTTGAGAAAAGT' #HCoV_229E
#trsL = 'TTTCGTTTAGTTGAGAA' #HCoV_NL63
#trsL = 'ATTTCGTTTAGTTCGA' #FCoV
#trsL = 'ATTTCGTTTAGTTCGA' #TGEV

#Threshold depends on mfe for canonical trs_l/b matches
threshold = #HCoV_229E
#threshold = #HCoV_NL63
#threshold = #FCoV
#threshold = #TGEV

sequence = ''    
header = ''
for record in SeqIO.parse(file, "fasta"):
    header = record.id
    sequence = str(record.seq).upper()


leaderMatches = {}
for window in sequence[k:k+len(trsL)]
    cmd = f"RNAduplex --noLP"
    #input needed: variables trsL and sequence
    subprocess.run(cmd.split(), stdout=outputStream, check=TRUE)
    
    #get mfe-value from RNAduplex result
    if mfe < threshold
        #add window substring to the list: 
        endPosition = k+len(window)       
        leaderMatches[f"match {k}-{endPosition}"] = window


#list entries written to file -> can then be used for mlocarna:
with open(outfile, 'w') as outputStream:
    for match in leaderMatches:
        matchSeq = leaderMatches[match]
        outputStream.write(f">{match}\n{matchSeq}\n")

