## Conserved RNA secondary structures in direct proximity of Coronavirus TRS-B sequences may contribute to discontinuous transcription

### Supplement Files and Scripts

Welcome!

This github repository is part of our publication [Lamkiewicz et al. (2019)](https://www.tobepublished.com) and contains
our Python scripts as well as a small example dataset for the Human Coronavirus 229E.

In order to get started, simply clone this repository on your local disk.

```bash
https://github.com/klamkiew/cov_trs_structure.git && cd cov_trs_structure/
```

#### Requirements

In order to run our Python scripts without any problems, please make sure that you have a working
[Python3](https://www.python.org/downloads/) installation on your system. Further, we developed
our scripts under a Linux-based environment, thus, we do not guarantee our code to run on
MacOS or Windows.

Further, we used some packages that are not built in the standard Python3. However, 
you can simply install these packages with the following command:

```bash
# make sure that the 'pip3' is actually connected to your Python3 installation.
pip3 install -r requirements.txt
```


We further use the function of `RNAcofold v 2.4.11` which is implemented in
the `ViennaRNA 2` package. Please follow this [link here](https://www.tbi.univie.ac.at/RNA/) 
and the instruction on the webpage. Make sure that the programs are visible in 
your `$PATH` variable. By default, this is done automatically when installing the
ViennaRNA package without any parameters. For local installations and other things,
we refer to the nice [guideline](https://www.tbi.univie.ac.at/RNA/#download) provided
by the authors of ViennaRNA.


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

`leader_hybridization.py` works in a similar fashion and contains the code to create the set of negative 
instances. Again, it expects the *reference* and *csv-file* as input parameters and further needs an *output file*
in order to store the extracted sequences.


```bash
python3 leader_hybridization.py example/HCoV_229E.fa example/hcov229e_trs.csv example/negative_seq_hcov229e.fa

```

#### Contact

If you have questions or comments, do not hesitate to contact us:

kevin.lamkiewicz{at}uni-jena{dot}de

#### License
GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

For more details, we refer to the [LICENSE](./LICENSE)
