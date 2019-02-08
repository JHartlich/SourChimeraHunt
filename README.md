# SourChimeraHunt

is used to extract chimera sequences from a multiple CCS containing FASTA file. 

## Usage
```bash
./SourChimeraHunt.py [INPUT FASTA] [DATABASE] [-Options]
```
### Positional Arguments:
```
  INPUT           name of input FASTA file containing multiple FASTA sequences
  DATABASE        name of database in FASTA format used for reference based chimera detection vie vsearch
```

### Optional Arguments:
```
  - h, --help     show help message and exit
  - sdb SECONDB   name of database in FASTA format used for secondary reference based chimera detection via vsearch
  - k   KSIZE     k-mer size for hash calculation, default: 20 nt
  - l   LENGTH    length cut off to filter to short sequences, default: 1400 nt
```
