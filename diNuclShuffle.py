#!/usr/bin/env python3

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

cmd = f"bash -c 'RNAalifold -r --cfactor 0.6 --nfactor 0.5 --noLP --noPS {ORIGINAL}'"
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=TRASH)
output = p.stdout.read().decode('ascii').split('\n')[1]
ENERGY = re.findall("(-\d+\.\d+)", output)[0]

energies = [ENERGY]

cmd = f"multiperm -n 1000 -w {ORIGINAL}"
p = subprocess.Popen(cmd, shell=True)
p.wait()

for data in glob("perm*.aln"):
  if os.path.exists(f"{DIR}/{data}"):
    os.remove(f"{DIR}/{data}")
  shutil.move(data, DIR)

for aln in glob(f"{DIR}/perm*aln"):
  cmd = f"bash -c 'RNAalifold -r --cfactor 0.6 --nfactor 0.5 --noLP --noPS {aln}'"
  p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=TRASH)
  output = p.stdout.read().decode('ascii').split('\n')[1]
  energies.extend(re.findall("(-\d+\.\d+)", output))
  os.remove(aln)
TRASH.close()  

a = np.array(energies).astype(float)
z = zscore(a)[0]
pvalue = norm.sf(abs(z)) * 2

print(f"MFE of original alignment: {ENERGY} kcal/mol")
print(f"p-value based on dinucl shuffle and z-score analysis: {pvalue}")
