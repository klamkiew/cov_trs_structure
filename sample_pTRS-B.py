#!/usr/bin/env python3

"""

TRS-L matching to genome to find putative hybridization subsequences.
This script is part of our publication

Lamkiewicz et al. (2021), "", Viruses, doi:



Usage:
    leader_hybridization.py [options] <GENOME_FASTA> <CSFILE> <OUTPUTFILE>

Options:
    -h, --help                          Show this little neat help message and exit.
    -r FLANKING, --region FLANKING      Size of extracted region up and downstream of the TRS_L. [default: 4]


Dependencies:
    BioPython 1.73 or greater
    docopt 0.6.2 or greater
    numpy 1.15.4 or greater

    ViennaRNA 2.4.13 installed in your $PATH variable.


Contact:
    kevin.lamkiewicz{at}uni-jena{dot}de
"""

###################################

import re
import subprocess
import sys

from Bio import SeqIO
import Levenshtein
from docopt import docopt
import numpy as np

###################################

d_complementNucl = { 
    "A":"T",
    "C":"G",
    "G":"C",
    "T":"A"
}

def reverseComplement(query):
    """
    Return the reverse complementary sequence of a given query.
    
    Parameters:
    query -- the given nucleotide sequence

    Return:
    String containing the reverse complementary nucleotide sequence

    """

    return(''.join([d_complementNucl[x] for x in query[::-1]]))

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

def apply_cofold(leader, fragment):
    """
    Wrapping function to call RNAcofold from the ViennaRNA package.
    Takes the two RNA molecules of interest and applies RNAcofold with the --noLP parameter.

    Parameters:
    leader -- first RNA molecule, the TRS-L sequence of our virus
    fragment -- second RNA molecule that interacts with the leader

    Return:
    minimum free energy value stored as float 

    """
    cmd = f'RNAcofold --noLP <<< "{leader}&{fragment}"'
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, executable="/bin/bash")
    cofoldResult = str(proc.stdout.read())

    # regex = re.compile(r'-?\d+.\d{2}')
    # the AttributeError should never occur, only if RNAcofold does not
    # terminate properly.
    try:
        mfe = float(re.search(regex,cofoldResult).group())
    except AttributeError:
        print(cofoldResult)
        exit(1)
    return mfe

###################################

if __name__ == "__main__":

    ###################################
    # Argument parsing

    args = docopt(__doc__)
    fileFasta = args['<GENOME_FASTA>']
    fileCS = args['<CSFILE>']
    outfile = args['<OUTPUTFILE>']
    flankingSize = int(args['--region'])

    ###################################
    # reading canonical core sequences

    s_canonicalCS = set()
    d_coreSequences = {}
    with open(fileCS, 'r') as inputStream:
        for line in inputStream:
            currentEntry = line.rstrip('\n').split(',')
            d_coreSequences[currentEntry[0]] = (currentEntry[1], int(currentEntry[2]))
            s_canonicalCS.add(int(currentEntry[2]))

    ###################################
    # reading the reference genome

    d_fastaRecords = {record.id : str(record.seq).upper() for record in SeqIO.parse(fileFasta,"fasta")}    
    header = ''
    sequence = ''
    if len(d_fastaRecords) != 1:
        print(f"Error: More than one fasta record found in {fileFasta}.")
        sys.exit(1)
    else:
        header = list(d_fastaRecords.keys())[0]
        sequence = d_fastaRecords[header]

    ###################################
    # extracting canonical interactions

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
        
    ###################################
    # calculate the interaction strength
    # of canonical energies.    

    canonicalEnergies = []
    regex = re.compile(r'-?\d+.\d{2}')
    for sgRNA, fragment in d_interactions.items():
        canonicalEnergies.append(apply_cofold(leader,fragment))
    
    ###################################
    # find regions in the genome
    # with similar interaction strength

    leaderCS = d_coreSequences['L'][0]
    canonicalEnergies = np.array(canonicalEnergies)
    canonicalTRS = [x[1] for x in d_coreSequences.values()]
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

    ###################################
    # write negative instances into
    # an fasta formatted file

    with open(outfile, 'w') as outputStream:
        for idx,sequence in enumerate(negativeSet):
            outputStream.write(f">negative_pseudoTRS_sequence_{idx+1}\n{sequence}\n")
            
