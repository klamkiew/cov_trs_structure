#!/usr/bin/env python3

"""

TRS sequence finder for Coronaviruses with flanking regions.

Usage:
    flanking_regions.py [options] <FASTA> <OUTPUTFILE>

Options:
    -h, --help                                  Show this little neat help message and exit.
    -r REGION, --region REGION                  Size of extracted region up- and downstream of the TRS. [default: 4]
"""

from Bio import SeqIO
import RNA
import sys
from docopt import docopt

def find_all(a_str, sub):
    return([k for k in range(len(a_str)) if a_str[k:k+len(sub)] == sub])    

args = docopt(__doc__)
#print(args)


file = args['<FASTA>']
outfile = args['<OUTPUTFILE>']
regionSize = int(args['--region'])

sequence = ''

trsSeqs = {
    "L": "TCTCAACT",
    "S": "TCTCAACT",
    "4": "ACTAAACT",
    "E": "TCTCAACT",
    "M": "TCTAAACT",
    "N": "TCTAAACT" 
}   

#trsSeqs = {
#    "L": "TCTCAACT"
#}

header = ''
for record in SeqIO.parse(file, "fasta"):
    header = record.id
    sequence = str(record.seq).upper()


trsRegion = { geneName : find_all(sequence, trs) for geneName, trs in trsSeqs.items() }
#for geneName, trs in trsSeqs.items():
#    trsRegion[geneName] = find_all(sequence, trs)

uniquePositions = set([x for liste in trsRegion.values() for x in liste])
#for liste in trsRegion.values():
#    for element in liste:
#        if not element in neueListe:
#            uniquePositions.append(element)

with open(outfile, 'w') as outputStream:
    for trs in uniquePositions:
        if trs - regionSize < 0:
            continue
        
        subSeq = sequence[trs-regionSize:trs+8+regionSize]
        outputStream.write(f">{header}|{trs-regionSize}-{trs+8+regionSize}\n{subSeq}\n")
        
        #print(f">{header}  {trs-80}-{trs}")
        #print(subSeq)
        #fc = RNA.fold_compound(subSeq)
        #structure, energy = fc.mfe()
        #print(fc.pf())
