## Conserved RNA secondary structures in direct proximity of Coronavirus TRS-B sequences may contribute to discontinuous transcription

### Supplement Files and Scripts

Welcome!

This github repository is part of our publication [Lamkiewicz et al. (2019)]{www.tobepublished.com} and contains
our Python scripts as well as a small example dataset for the Human Coronavirus 229E.

#### Python Scripts

In `find_trs.py` you'll find the code to extract the core sequences and their corresponding
upstream regions within the genome. You need a *reference genome* and a *csv-file* containing
positional information of the canonical TRS sequences.
In our example, the information for HCoV-229E is stored in `example/hcov229e_trs.csv`:

```
L,TCTCAACT,61
S,TCTCAACT,20554
4,ACTAAACT,24046
E,TCTCAACT,24587
M,TCTAAACT,24978
N,TCTAAACT,25667
```
The first column contains the *name* of the gene, the second column in turn the *core sequence*.
The third column is not needed for `find_trs.py`, but it contains the *starting position* of the core sequence
within the reference genome.

In order to extract upstream regions, you simply want to call `find_trs.py` with the reference, the csv-file
and some *output file* which will store the extracted sequences.

```bash
# if -t is activated, the upstream region plus the canonical TRS sequence is extracted
python3 find_trs.py [-t] example/HCoV_229E.fa example/hcov229e_trs.csv example/extracted_seq_hcov229e.fa

```

