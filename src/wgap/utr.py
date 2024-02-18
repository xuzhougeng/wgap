import pyranges as pr
from typing import List
from .model_loader import gff3_loader
from pyranges import PyRanges



def add_utr(gff_file: str, gtf_file_list: List[str]):
    gene_dict = gff3_loader(gff_file)