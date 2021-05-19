#!/usr/bin/env python3

"""

Wrapper for multiperm [1] and RNAalifold [2] to check for significance of an
consensus RNA structure of interest.

The input alignment in CLUSTAL format is folded with RNAalifold
to get the original MFE. Then, 1000 shuffled alignments are created
with multiperm (approximate dinucleotide shuffling of alignments).
For each of these alignments RNAalifold is invoked again.
A simple z-score analysis and p-value calculation is then applied
to check whether the original consensus structure is significant
considering the energy with (di)nucleotide context.

Usage:
  diNuclShuffle.py <CLUSTALW>

Dependencies:
    numpy>=1.15.4
    scipy>=1.6.3


    ViennaRNA 2.4.13 installed in your $PATH variable.
    multiperm installed in your $PATH variable.

Contact:
  kevin.lamkiewicz{at}uni-jena{dot}de

"""

import os
import shutil
import subprocess
import re
import sys

from glob import glob

import numpy as np
from scipy.stats import norm
from scipy.stats.mstats import zscore


ORIGINAL = sys.argv[1]
DIR = os.path.dirname(ORIGINAL)

TRASH = open(os.devnull)


###################################
# original MFE calculation
cmd = f"bash -c 'RNAalifold -r --cfactor 0.4 --nfactor 0.5 --noLP --noPS {ORIGINAL}'"
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=TRASH)
output = p.stdout.read().decode('ascii').split('\n')[1]
ENERGY = re.findall("(-\d+\.\d+)", output)[0]

energies = [ENERGY]

###################################
# Every day I'm shufflin'
cmd = f"multiperm -n 1000 -w {ORIGINAL}"
p = subprocess.Popen(cmd, shell=True)
p.wait()

###################################
# move tmp alignments to the input directory
for data in glob("perm*.aln"):
  if os.path.exists(f"{DIR}/{data}"):
    os.remove(f"{DIR}/{data}")
  shutil.move(data, DIR)

###################################
# for each shuffled alignment
# determine the MFE value

for aln in glob(f"{DIR}/perm*aln"):
  cmd = f"bash -c 'RNAalifold -r --cfactor 0.6 --nfactor 0.5 --noLP --noPS {aln}'"
  p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=TRASH)
  output = p.stdout.read().decode('ascii').split('\n')[1]
  energies.extend(re.findall("(-\d+\.\d+)", output))
  os.remove(aln)  # removing all tmp alignments.
TRASH.close()

###################################
# calculate z-score and p-values
a = np.array(energies).astype(float)
z = zscore(a)[0]
pvalue = norm.sf(abs(z)) * 2

print(f"MFE of original alignment: {ENERGY} kcal/mol")
print(f"p-value based on dinucl shuffle and z-score analysis: {pvalue}")
