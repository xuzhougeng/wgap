# WGAP: Whole Genome Annotation Pipeline


## Installation

You can install WGAP using bioconda

```bash
conda create -n wgap -c bioconda wgap
```

Alternatively, you can also install the dependency manually and then install the wgap from github

```bash
git clone 
cd 
python setup.py install
```

## Usage

WGAP suppory a widely 

- transcripts.fasta
- transcripts.gff
- protein.fa
- protein.gff
- snap_hmm
- augustus_species

### Example 1: run with fasta as evidence and exsiting gene models

Firstly, generate the configuration file

```bash
wgap init --transcripts_fasta transcripts.fasta --protein_fasta protein.fa
```




## 开发计划