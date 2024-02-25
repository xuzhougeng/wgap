# download the protin from UniPort

import gzip
import urllib.request
from Bio import SwissProt


def download_uniprot(specie, dataset):
    #wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_plants.dat.gz
    #zcat uniprot_sprot_plants.dat.gz |\
    #    awk '{if (/^ /) {gsub(/ /, ""); print} else if (/^AC/) print ">" $2}' |\
    #    sed 's/;$//'> protein.fa
    file_name = "uniprot_{dataset}_{specie}.dat.gz".format(dataset=dataset, specie=specie)
    file_path = file_name
    base_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/"
    url = base_url + file_name

    try:
        urllib.request.urlretrieve(url, filename=file_path)
    except:
        print("You can also download protein with `wget {}`, and then rerun wgap".format(url))
    
    return file_path


# convert SwissPort to FASTA
def convet_dat_to_fasta(dat, fasta):
    handle = gzip.open(dat, "rt")
    file_out = open(fasta, "w")
    for record in SwissProt.parse(handle):
        gene_name = record.gene_name
        organism = record.organism
        seq = record.sequence
        out_line = ">{}{}\n{}\n".format(gene_name, organism, seq)
        file_out.writelines(out_line)


def add_exon_to_gene_model(input_file, output_file):
    # 读取文件内容
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # 处理每一行，准备修改后的内容
    modified_lines = []
    for line in lines:
        # 跳过注释行
        if line.startswith('#'):
            modified_lines.append(line)
            continue

        # 写入原始行
        modified_lines.append(line)

        # 处理 CDS 行
        if '\tCDS\t' in line:
            fields = line.strip().split('\t')
            fields[2] = 'exon'  # 替换类型为 exon
            fields[8] = fields[8].replace('CDS', 'exon')  # 替换 ID 中的 CDS 为 exon
            new_line = '\t'.join(fields) + '\n'
            modified_lines.append(new_line)

    # 写回同一个文件
    with open(output_file, 'w') as outfile:
        outfile.writelines(modified_lines)

    

def fix_fasta(fasta, out_file):
    # MAKER only support A, C, G, T, N
    # This script will convert all other characters(RYSWKMBDHV) to N
    base = ['A', 'T', 'C', 'G', 'N']

    with open(out_file, 'w' ) as handler:
        for line in open(fasta, 'r') :
            if line.startswith('>'):
                handler.write(line)
                continue

            new_line = ''
            for char in line.strip():
                if char.upper() not in  base:
                    new_line += 'N'
                else:
                    new_line += char
            handler.write(new_line + '\n')

def parse_gff3(gff_file):
    genes = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            feature_type, strand = parts[2], parts[6]  # 获取特征类型和链方向
            if feature_type not in ['exon', 'CDS']:
                continue

            attributes = parts[8]
            gene_id = None
            for attr in attributes.split(';'):
                if attr.startswith('Parent='):
                    gene_id = attr.split('=')[1]
                    break
            if not gene_id:
                continue

            start, end = int(parts[3]), int(parts[4])

            if gene_id not in genes:
                genes[gene_id] = {'exon': [], 'CDS': [], 'strand': strand}  # 记录链方向
            genes[gene_id][feature_type].append((start, end))

    return genes

def calculate_utr_lengths(genes):
    utr_lengths = {}
    for gene_id, features in genes.items():
        exons = sorted(features['exon'], key=lambda x: x[0])
        CDS = sorted(features['CDS'], key=lambda x: x[0])
        strand = features['strand']  # 获取链方向

        if not CDS:  # 如果没有CDS信息，跳过
            continue

        # 根据链方向调整5'UTR和3'UTR的计算
        if strand == '+':
            five_prime_utr_length = calculate_five_prime_utr_length(exons, CDS[0][0])
            three_prime_utr_length = calculate_three_prime_utr_length(exons, CDS[-1][1])
        else:
            five_prime_utr_length = calculate_three_prime_utr_length(exons, CDS[-1][1], reverse=True)
            three_prime_utr_length = calculate_five_prime_utr_length(exons, CDS[0][0], reverse=True)

        utr_lengths[gene_id] = (five_prime_utr_length, three_prime_utr_length)

    return utr_lengths

def calculate_five_prime_utr_length(exons, cds_start, reverse=False):
    utr_length = 0
    for exon_start, exon_end in (reversed(exons) if reverse else exons):
        if exon_end < cds_start:
            utr_length += exon_end - exon_start + 1
        elif exon_start < cds_start:
            utr_length += cds_start - exon_start
            break
    return utr_length

def calculate_three_prime_utr_length(exons, cds_end, reverse=False):
    utr_length = 0
    for exon_start, exon_end in (exons if reverse else reversed(exons)):
        if exon_start > cds_end:
            utr_length += exon_end - exon_start + 1
        elif exon_end > cds_end:
            utr_length += exon_end - cds_end
            break
    return utr_length

def utr_stat(gff_file,output_file):
    genes = parse_gff3(gff_file)
    utr_lengths = calculate_utr_lengths(genes)

    # 总基因数量
    total_genes = len(genes)

    # 具有UTR的基因数量
    genes_with_utr = sum(1 for lengths in utr_lengths.values() if lengths[0] > 0 or lengths[1] > 0)

    # 计算总的5'UTR和3'UTR长度
    total_five_prime_utr_length = sum(lengths[0] for lengths in utr_lengths.values())
    total_three_prime_utr_length = sum(lengths[1] for lengths in utr_lengths.values())

    # 计算平均UTR长度，如果有具有UTR的基因
    average_five_prime_utr_length = total_five_prime_utr_length / genes_with_utr if genes_with_utr > 0 else 0
    average_three_prime_utr_length = total_three_prime_utr_length / genes_with_utr if genes_with_utr > 0 else 0

    with open(output_file, 'w') as f:
        f.write(f"Total genes: {total_genes}\n")
        f.write(f"Genes with UTR: {genes_with_utr} ({genes_with_utr / total_genes * 100:.2f}%)\n")
        f.write(f"Average 5'UTR length: {average_five_prime_utr_length:.2f}\n")
        f.write(f"Average 3'UTR length: {average_three_prime_utr_length:.2f}\n\n")

        for gene_id, lengths in utr_lengths.items():
            f.write(f"{gene_id}\t5'UTR: {lengths[0]}\t3'UTR: {lengths[1]}\n")