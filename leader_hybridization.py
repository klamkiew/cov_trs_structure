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
import numpy as np
import Levenshtein

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

def apply_cofold(leader, fragment):
    cmd = f'RNAcofold --noLP <<< "{leader}&{fragment}"'
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, executable="/bin/bash")
    cofoldResult = str(proc.stdout.read())
    try:
        mfe = float(re.search(regex,cofoldResult).group())
    except AttributeError:
        print(cofoldResult)
        exit(1)
    return mfe

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
    regex = re.compile(r'-?\d+.\d{2}')
    for sgRNA, fragment in d_interactions.items():
        canonicalEnergies.append(apply_cofold(leader,fragment))
        print(sgRNA, d_coreSequences[sgRNA], fragment, canonicalEnergies[-1])
    
    leaderCS = d_coreSequences['L'][0]
    canonicalEnergies = np.array(canonicalEnergies)
    canonicalTRS = [x[1] for x in d_coreSequences.values()]
    #print(canonicalTRS)
    ranges = [range(x-50,x+50) for x in canonicalTRS]
    
    negativeSamplingStart = d_coreSequences['L'][1] + len(leaderCS) + flankingSize
    negativeSet = []
    for i in range(negativeSamplingStart, len(sequence)-len(leaderCS)):
        fragment = reverseComplement(sequence[i:i+len(leader)])

        mfe = apply_cofold(leader,fragment)
        if mfe < (np.mean(canonicalEnergies) + np.std(canonicalEnergies)):
            negativeSequence = sequence[i-150: i+len(leaderCS)]
            if any([Levenshtein.ratio(negativeSequence, x) > 0.9 for x in negativeSet]) or any([i in x for x in ranges]):
                continue
            negativeSet.append(negativeSequence)
            print(f"pTRS-B\_{len(negativeSet)} & {i-flankingSize}--{i+len(leaderCS)+flankingSize} & {leader.replace('T','U')} & {fragment.replace('T','U')} & {mfe}\,kcal/mol \\\\")
            print(f"{i-flankingSize}--{i+len(leaderCS)+flankingSize}")
    with open(outfile, 'w') as outputStream:
        #constraint = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        for idx,sequence in enumerate(negativeSet):
            outputStream.write(f">negative_pseudoTRS_sequence_{idx+1}\n{sequence}\n")
            #outputStream.write(f"{constraint[:len(leaderCS)].rjust(len(sequence),'.')} #1\n")
