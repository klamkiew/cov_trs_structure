#!/usr/bin/env python3

"""

TRS sequence finder for Coronaviruses!

Usage:
    find_trs.py [options] <FASTA> <OUTPUTFILE>

Options:
    -h, --help                                  Show this little neat help message and exit.
    -r REGION, --region REGION                  Size of extracted region downstream of the TRS. [default: 100]
    -m MISMATCH, --mismatch MISMATCH            Allow one mismatch in trs sequence if set to 1. [default: 0]
    -trs INCLUDE_TRS, --withtrs INCLUDE_TRS     Include the trs sequence in extracted sequence if set to 1. [default: 0]
"""

from Bio import SeqIO
import RNA
import sys
from docopt import docopt

def find_all(a_str, sub):
    return([k for k in range(len(a_str)) if a_str[k:k+len(sub)] == sub])    

def hamming(s1, s2):
    if len(s1) == len(s2):
        return int(sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2)))

def find_with_mism(a_str, sub):   #find match with one mismatch allowed
    #return([k for k in range(len(a_str)) if hamming(a_str[k:k+len(sub)], sub) == 1 or hamming(a_str[k:k+len(sub)], sub) == 0]) 

    indexes = []
    k = 0
    while k < len(a_str):
        if hamming(a_str[k:k+len(sub)], sub) == 1 or hamming(a_str[k:k+len(sub)], sub) == 0:  
            indexes.append(k)
        k += 1
    return indexes

args = docopt(__doc__)
#print(args)


file = args['<FASTA>']
outfile = args['<OUTPUTFILE>']
regionSize = int(args['--region'])
mismatch = int(args['--mismatch'])
withtrs = int(args['--withtrs'])
#if not len(sys.argv) > 3:
#    regionSize = 100
#else: 
#    regionSize = sys.argv[3]


sequence = ''

trsSeqs = {
    "S": "TCTCAACT",
    "4": "ACTAAACT",
    "E": "TCTCAACT",
    "M": "TCTAAACT",
    "N": "TCTAAACT"
}   

header = ''
for record in SeqIO.parse(file, "fasta"):
    header = record.id
    sequence = str(record.seq)

if mismatch == 0:
    trsRegion = { geneName : find_all(sequence, trsB) for geneName, trsB in trsSeqs.items() }
    print(trsRegion)
    #for geneName, trsB in trsSeqs.items():
    #    trsRegion[geneName] = find_all(sequence, trsB)
elif mismatch == 1:
    trsRegion = { geneName : find_with_mism(sequence, trsB) for geneName, trsB in trsSeqs.items() }

uniquePositions = set([x for liste in trsRegion.values() for x in liste])
#for liste in trsRegion.values():
#    for element in liste:
#        if not element in neueListe:
#            uniquePositions.append(element)

with open(outfile, 'w') as outputStream:
    for trs in uniquePositions:
        if trs - regionSize < 0:
            continue
        
        if withtrs == 0:
            subSeq = sequence[trs-regionSize:trs]
            outputStream.write(f">{header}|{trs-regionSize}-{trs}\n{subSeq}\n")
        
        #option for subsequence including trs sequence:
        if withtrs == 1:
            subSeq = sequence[trs-regionSize:trs+8]
            outputStream.write(f">{header}|{trs-regionSize}-{trs+8}\n{subSeq}\n")
        
        #print(f">{header}  {trs-80}-{trs}")
        #print(subSeq)
        #fc = RNA.fold_compound(subSeq)
        #structure, energy = fc.mfe()
        #print(fc.pf())


