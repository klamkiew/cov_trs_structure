#!/usr/bin/env python3

"""

TRS-L matching to genome to find putative hybridization subsequences.

Usage:
    leader_hybridization.py [options] <GENOME_FASTA> <CSFILE> <OUTPUTFILE>

Options:
    -h, --help                          Show this little neat help message and exit.
    -r FLANKING, --region FLANKING      Size of extracted region up and downstream of the TRS_L. [default: 4]

"""

from Bio import SeqIO
import RNA
import sys
import csv
import subprocess
from docopt import docopt
import re

d_complementNucl = { 
    "A":"T",
    "C":"G",
    "G":"C",
    "T":"A"
}

def reverseComplement(query):
    return(''.join([d_complementNucl[x] for x in query[::-1]]))

def find_all(a_str, sub):
    return([k for k in range(len(a_str)) if a_str[k:k+len(sub)] == sub])    


if __name__ == "__main__":
    
    args = docopt(__doc__)
    fileFasta = args['<GENOME_FASTA>']
    fileCS = args['<CSFILE>']
    outfile = args['<OUTPUTFILE>']
    flankingSize = int(args['--region'])
    s_canonicalCS = set()
    d_coreSequences = {}
    with open(fileCS, 'r') as inputStream:
        for line in inputStream:
            currentEntry = line.rstrip('\n').split(',')
            d_coreSequences[currentEntry[0]] = (currentEntry[1], int(currentEntry[2]))
            s_canonicalCS.add(int(currentEntry[2]))


    d_fastaRecords = {record.id : str(record.seq).upper() for record in SeqIO.parse(fileFasta,"fasta")}    
    header = ''
    sequence = ''
    if len(d_fastaRecords) != 1:
        print(f"Error: More than one fasta record found in {fileFasta}.")
        sys.exit(1)
    else:
        header = list(d_fastaRecords.keys())[0]
        sequence = d_fastaRecords[header]

    d_interactions = {}
    leader = ''
    for sgRNA, (coreSequence, position) in d_coreSequences.items():
        occurences = find_all(sequence, coreSequence)
        
        if position not in occurences:
            print(f"Error! Couldn't find {coreSequence} of sg mRNA {sgRNA} in the reference.")
            sys.exit(2)
        
        fragment = sequence[position - flankingSize : position + len(coreSequence) + flankingSize]
        if sgRNA == "L":
            leader = fragment
        else:
            d_interactions[sgRNA] = reverseComplement(fragment)
        
    canonicalEnergies = []
    regex = re.compile(r'-\d+.\d{2}')
    for sgRNA, fragment in d_interactions.items():
        cmd = f'RNAcofold --noLP <<< "{leader}&{fragment}"'
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, executable="/bin/bash")
        cofoldResult = str(proc.stdout.read())
        canonicalEnergies.append(float(re.search(regex,cofoldResult).group()))
    
    negativeSamplingStart = #leader pos + cs length
    for i in range()
#trsL = 'CGTTTAGTTGAGAAAAGT' #HCoV_229E
#trsL = 'TTTCGTTTAGTTGAGAA' #HCoV_NL63
#trsL = 'ATTTCGTTTAGTTCGA' #FCoV
#trsL = 'ATTTCGTTTAGTTCGA' #TGEV

#Threshold depends on mfe for canonical trs_l/b matches
#threshold = #HCoV_229E
#threshold = #HCoV_NL63
#threshold = #FCoV
#threshold = #TGEV

# sequence = ''    
# header = ''
# for record in SeqIO.parse(file, "fasta"):
#     header = record.id
#     sequence = str(record.seq).upper()


# leaderMatches = {}
# for window in sequence[k:k+len(trsL)]:
#     cmd = f"RNAduplex --noLP"
#     #input needed: variables trsL and sequence
#     subprocess.run(cmd.split(), stdout=outputStream, check=TRUE)
    
#     #get mfe-value from RNAduplex result
#     if mfe < threshold:
#         #add window substring to the list: 
#         endPosition = k+len(window)       
#         leaderMatches[f"match {k}-{endPosition}"] = window


# #list entries written to file -> can then be used for mlocarna:
# with open(outfile, 'w') as outputStream:
#     for match in leaderMatches:
#         matchSeq = leaderMatches[match]
#         outputStream.write(f">{match}\n{matchSeq}\n")

