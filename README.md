# SourChimeraHunt
is a command line based tool to extract chimera sequences from a multiple CCS containing FASTA file.\
It uses MinHashing from [sourmash](https://github.com/dib-lab/sourmash "sourmash @ Github") Python API ([Brown et al. 2016](https://joss.theoj.org/papers/3d793c6e7db683bee7c03377a4a7f3c9)) and reference based chimera detection from [vsearch](https://github.com/torognes/vsearch "vsearch @ GitHub") ([Rognes et al. 2016](https://peerj.com/articles/2584/)) to do so.\
To operate SourChimeraHunt needs at least the positional arguments: the input file and a reference database, both in FASTA format.


## Requirements
[Python](https://www.python.org/downloads "Download Python") 2.7 or 3\
[sourmash](https://github.com/dib-lab/sourmash "sourmash @ Github")\
[vsearch](https://github.com/torognes/vsearch "vsearch @ GitHub")


## Usage
```bash
./SourChimeraHunt.py [INPUT FASTA] [DATABASE] [-Options]
```
### Positional Arguments:
```
  INPUT           name of input FASTA file containing multiple FASTA sequences
  DATABASE        name of database in FASTA format used for reference based chimera detection via vsearch
```

### Optional Arguments:
```
  - h, --help     show help message and exit
  - sdb SECONDB   name of database in FASTA format used for secondary reference based chimera detection via vsearch
  - k   KSIZE     k-mer size for hash calculation, default: 20 nt
  - l   LENGTH    length cut off to filter to short sequences, default: 1400 nt
```
