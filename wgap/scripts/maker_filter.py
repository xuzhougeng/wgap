#!/usr/bin/env python3

"""
filter the maker output gff by:
    - AED
    - eAED
    - QI
    - mRNA length and mRNA number
"""

from math import inf
import numpy as np
import re
import argparse


def parse_anno(col):
    """
    parse the annotation column, and return a dict
    """
    anno = re.split('[;=]', col)
    if anno[-1] == '':
        anno.pop()

    anno_dict = {}
    for i in range(0, len(anno), 2):
        anno_dict[anno[i]] = anno[i+1]
    return anno_dict



class QI:
    """
    MAKER QI value
    """
    def __init__(self, exon_num, prt_len,
                five_utr_len, three_utr_len,
                fct_ss_conf_est, fct_exon_ovlp_est, fct_exon_ovlp_prt,
                fct_ss_conf_ab, fct_exon_ovlp_ab):

        self.exon_num = int(exon_num)                     # Number of exons in the mRNA
        self.prt_len = int(prt_len)                       # Length of the protein sequence produced by the mRNA (-l)
        self.five_utr_len  = int(five_utr_len)            # Length of the 5' UTR,
        self.three_utr_len = int(three_utr_len)           # Length of the 3' UTR
        self.fct_ss_conf_est = float(fct_ss_conf_est)     # Fraction of splice sites confirmed by an EST alignment
        self.fct_exon_ovlp_est = float(fct_exon_ovlp_est) # Fraction of exons that overlap an EST alignmetn(-e)
        self.fct_exon_ovlp_prt = float(fct_exon_ovlp_prt) # Fraction of exons that overlap EST or Protein alignments(-o)
        self.fct_ss_conf_ab = float(fct_ss_conf_ab)       # Fraction of splice site confrimed by a ab-initio prediction(-a)
        self.fct_exon_ovlp_ab = float(fct_exon_ovlp_ab)   # Fraction of exons that overlap a ab-initio prediction(-t)
    
    @classmethod
    def parse_qi(cls, line):
        """parse the MAKER QI record
        """
        #_QI=170|1|1|1|0|0|3|183|340
        items = line.strip().split('|')
        
        five_utr_len = items[0]
        fct_ss_conf_est = items[1]
        fct_exon_ovlp_est = items[2]
        fct_exon_ovlp_prt = items[3]
        fct_ss_conf_ab = items[4]
        fct_exon_ovlp_ab = items[5]
        exon_num = items[6]
        three_utr_len = items[7]
        prt_len = items[8]

        return cls( exon_num, prt_len,
                five_utr_len, three_utr_len,
                fct_ss_conf_est, fct_exon_ovlp_est, fct_exon_ovlp_prt,
                fct_ss_conf_ab, fct_exon_ovlp_ab)

    @staticmethod
    def comparison(qry, sub):
        """
        :param qry : QI to compare
        :param sub : QI to be compared
        """
        if qry.exon_num < sub.exon_num:
            return False
        if qry.prt_len < sub.prt_len:
            return False
        if qry.five_utr_len < sub.five_utr_len:
            return False
        if qry.three_utr_len < sub.three_utr_len:
            return False
        if qry.fct_ss_conf_est < sub.fct_ss_conf_est:
            return False
        if qry.fct_exon_ovlp_est < sub.fct_exon_ovlp_est:
            return False
        if qry.fct_exon_ovlp_prt < sub.fct_exon_ovlp_prt:
            return False
        if qry.fct_ss_conf_ab < sub.fct_ss_conf_ab:
            return False
        if qry.fct_exon_ovlp_ab < sub.fct_exon_ovlp_ab:
            return False
        
        return True

class mRNAInfo:
    
    def __init__(self, ID, len, parent=None, QI=None,AED=None):
        self.ID = ID
        self.len = len # from mRNA start to mRNA end with introns
        
        self.parent = parent
        self.QI = QI
        self.AED = AED

    def __str__(self):
        return self.ID
    
    def __repr__(self) -> str:
        return self.ID


class GeneInfo:

    def __init__(self, ID, len):
        self.ID = ID
        self.len = len
        self.mRNA_dict = {}

    def __str__(self) -> str:
        return self.ID
    
    def __repr__(self) -> str:
        return self.ID

    def add_mRNA(self, mRNA: mRNAInfo):
        self.mRNA_dict[mRNA.ID] = mRNA

def parse_gff(gff):
    """parse the gff file, and extract the gene model information
    :param gff: the gff file
    """
    gene_dict = {}

    for line in open(gff, "r"):
        if line.startswith("#"):
            continue

        line = line.strip()
        cols = line.split('\t')
        if len(cols) != 9:
            continue
        if cols[1] != 'maker':
            continue

        if cols[2] == 'gene':
            anno = parse_anno(cols[8])
            gene_id = anno['ID']
            gene_len = int(cols[4]) - int(cols[3])
            gene_dict[gene_id] = GeneInfo(gene_id, gene_len)
        elif cols[2] == 'mRNA':
            anno = parse_anno(cols[8])
            QI_ = anno['_QI']
            mRNA_ID = anno['ID']
            mRNA_len = int(cols[4]) - int(cols[3])

            mRNA = mRNAInfo(mRNA_ID, mRNA_len)
            mRNA.QI = QI.parse_qi(QI_)
            mRNA.parent = anno['Parent']
            mRNA.AED = float(anno['_AED'])

            gene_dict[mRNA.parent].add_mRNA(mRNA)
        
    return gene_dict


def get_gene_len_list(gene_dict):
    gene_len_list = []
    for gene in gene_dict.values():
        gene_len_list.append(gene.len)
    return gene_len_list

def get_exon_num_list(gene_dict):
    exon_num_list = []
    for gene in gene_dict.values():
        for mRNA in gene.mRNA_dict.values():
            exon_num_list.append(mRNA.QI.exon_num)
    return exon_num_list

def filter_gff(gff, gene_dict, AED_thresh, QI_thresh, geneid, filter_outlier=False):
    """filter the gff by the gene list and mRNA list
    :param gff : the input gff file
    :param gene_dict : gene dictory
    :param AED_thresh : the AED threshhold
    :param QI_thresh : the QI threshhold
    :param geneid : the gene id list
    """

    min_exon_num = 0
    max_exon_num = inf
    
    min_gene_len = 0
    max_gene_len = inf
    if filter_outlier:
        gene_len_list = get_gene_len_list(gene_dict)
        exon_num_list = get_exon_num_list(gene_dict)

        min_exon_num,max_exon_num = np.percentile(exon_num_list, [25,75])
        min_gene_len,max_gene_len = np.percentile(gene_len_list, [25,75])

    good_gene,good_mRNA = set(),set()
    for gene in gene_dict.values():
        
        if gene.ID in geneid or gene.len < min_gene_len or gene.len > max_gene_len:
            continue

        for mRNA in gene.mRNA_dict.values():
            if mRNA.AED < AED_thresh or mRNA.QI.exon_num < min_exon_num or mRNA.QI.exon_num > max_exon_num:
                continue

            if QI.comparison(mRNA.QI, QI_thresh):
                good_mRNA.add(mRNA.ID)
                good_gene.add(gene.ID)

    for line in open(gff, "r"):
        if line.startswith("#"):
            continue
        line = line.strip()
        cols = line.split('\t')
        if len(cols) != 9:
            continue
        if cols[1] != 'maker':
            continue

        anno_dict = parse_anno(cols[8])

        # filter gene row by good_gene, filter exon/mRNA row by good_mRNA
        if cols[2] == 'gene':
            if anno_dict['ID'] in good_gene:
                print(line)
            else:
                continue
        elif cols[2] == 'mRNA':
            if anno_dict['ID'] in good_mRNA:
                print(line)
            else:
                continue
        else:
            if anno_dict['Parent'] in good_mRNA:
                print(line)
            else:
                continue


def get_opt():
    group = argparse.ArgumentParser()
    group.add_argument("-c", "--ss", type = float, default = -1, 
            help = "The faction of splice sites confirmed by an EST alignment, default -1", required = False)
    group.add_argument("-e", "--exon", type = float, default = -1, 
            help = "The faction of exons that overlap an EST alignment, default -1", required = False)
    group.add_argument("-o", "--exon2", type = float, default = -1, 
            help = "The faction of exons that overlap an EST/Protein alignment, default -1", required = False)
    group.add_argument("-a", "--ss2", type = float, default = -1, 
            help = "The faction of splice sites confirmed by an ab-initio prediction, default -1", required = False)
    group.add_argument("-t", "--exon3", type = float, default = -1, 
            help = "The faction of exon confirmed by an ab-initio prediction, default -1", required = False)
    group.add_argument("-l", "--length", type = int, default = 0, 
            help = "The min length of the protein sequence produced by the mRNA", required = False)
    group.add_argument("-d", "--AED", type = float, default = 1, required = False,
            help = "Max AED  to allow, default is 1")
    group.add_argument("-i", "--geneid", required = False,
            help = "filter by given gene list")
    group.add_argument("--filter-outlier", action = "store_true", default = False)
    group.add_argument("gff", help = "input gff file")
    
    return group.parse_args()

if __name__ == "__main__":
    opts = get_opt()
    gff = opts.gff

    AED_thresh = opts.AED
    QI_thresh = QI(0, opts.length, 0, 0, opts.ss, opts.exon, opts.exon2, opts.ss2, opts.ss2)

    gene_dict = parse_gff(gff)

    if opts.geneid is None:
        geneid = []

    filter_gff(gff, gene_dict, AED_thresh, QI_thresh, geneid, opts.filter_outlier)

