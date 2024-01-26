import pyranges as pr
from pyranges.pyranges_main import PyRanges
import os
import shutil

from typing import Dict, List, Union
from wgap.scripts.seq_utils import read_fasta, write_fasta, SeqRecordDict

def read_repeat_masker(repeat_file: str, skip :int = 3) -> PyRanges:
    
    # open file and skip header
    repeat_file_handle = open(repeat_file, 'r')
    for i in range(skip):
        next(repeat_file_handle)

    chromosome_list = []
    start_list = []
    end_list = []
    

    for line in repeat_file_handle:
        # Skip header and empty lines
        line = line.strip()
        if line.startswith('SW')  or line == '':
            continue

        # Extract fields from the line
        fields = line.split()

        # Parse fields and create Repeat object 
        chrom = fields[4]
        start = int(fields[5]) - 1
        end = int(fields[6])

        chromosome_list.append(chrom)
        start_list.append(start)
        end_list.append(end)
    
    return  pr.PyRanges(chromosomes = chromosome_list, starts = start_list, ends= end_list)

def split_genome(fasta: SeqRecordDict | str, repeat_masker: str, stringtie: str, homology:str=None,  out_dir:str="./") -> None:
    if isinstance(fasta, str):
        fasta = read_fasta(fasta)
    # fasta_pr
    fasta_pr = pr.PyRanges(chromosomes = list(fasta.keys()), starts = [0]*len(fasta), ends = [len(fasta[chrom]) for chrom in fasta])

    te_pr = read_repeat_masker(repeat_masker)
    # read stringtie
    stringtie_pr = pr.read_gtf(stringtie)
    stringtie_pr = stringtie_pr[stringtie_pr.Feature=="transcript"]

    # find overlap and filter
    te_pr = te_pr.overlap(stringtie_pr, invert=True, strandedness=False)
   
    # read homology
    if homology is not None:
        homology_pr = pr.read_gff3(homology)
        homology_pr = homology_pr[homology_pr.Feature=="mRNA"]
        te_pr = te_pr.overlap(homology_pr, invert=True, strandedness=False)
    
    # subtract the overlap
    fasta_pr2 = fasta_pr.subtract(te_pr)
    # remove the small interval less than 300
    fasta_pr2.Length = fasta_pr2.End - fasta_pr2.Start
    fasta_pr3 = fasta_pr2[fasta_pr2.Length > 300]

    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    os.makedirs(out_dir)
    
    
    for (chrom,start,end) in zip(fasta_pr3.Chromosome, fasta_pr3.Start, fasta_pr3.End):
        if (end -start ) < 300 :
            continue
        
        out_file = f"{out_dir}/{chrom}_{start}_{end}.fa"
        chrom_name = f"{chrom}_{start}_{end}"
        sequence = fasta[chrom][start:end]
        write_fasta({chrom_name: sequence}, out_file)

def split_genome2(fasta: SeqRecordDict | str, windows_size:int = 50000, out_dir:str="./") -> None:
    
    if isinstance(fasta, str):
        fasta = read_fasta(fasta)

    for chrom, sequence in fasta.items():
        for i in range(0, len(sequence), windows_size):
            start = i
            end = i + windows_size
            if end > len(sequence):
                end = len(sequence)
            out_file = f"{out_dir}/{chrom}_{start}_{end}.fa"
            chrom_name = f"{chrom}_{start}_{end}"
            write_fasta({chrom_name: sequence[start:end]}, out_file)

if __name__ == "__main__":
    import argparse
    args = argparse.ArgumentParser(description="Split genome into windows")
    args.add_argument("fasta", help="fasta file")
    # repeat masker result
    args.add_argument("-rm", "--repeat_masker", help="repeat masker result")
    # homology result
    args.add_argument("-hm", "--homology", help="homology result gff by miniprot")
    # stingtie result
    args.add_argument("-st", "--stingtie", help="stingtie result gff")
    # outdir
    args.add_argument("-o", "--outdir", help="output dir", default="./")

    args = args.parse_args()

    split_genome(args.fasta, args.repeat_masker, args.stingtie, args.homology, args.outdir)