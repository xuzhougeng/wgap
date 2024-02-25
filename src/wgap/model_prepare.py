from typing import Optional, Dict, List
from .seq_utils import SeqRecordDict, read_fasta
from .gene import Gene
from .model_loader import stringtie_loader
import numpy as np

from subprocess import Popen, PIPE

def calculate_mean_and_sigma(numbers):
    # 计算平均值和标准差
    mean = np.mean(numbers)
    sigma = np.std(numbers)
    return mean, sigma

def diamond_blastp(query, subject, threads=96, outfmt=6, out=None):
    """
    Run diamond blastp
    """
    # build diamond database
    diamond_db_cmd = f"diamond makedb --in {subject} -d {subject}"
    p = Popen(diamond_db_cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    # diamond_cmd = f"diamond blastp -q {query} -d {subject}.dmnd -p {threads} -f {outfmt} -o {out}"
    diamond_cmd = f"diamond blastp -q {query} -d {subject}.dmnd -p {threads} -f {outfmt}"
    p = Popen(diamond_cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    return stdout, stderr

def ncbi_blastp(query, subject, threads=96, outfmt=6, out=None):
    """
    Run blastp
    """
    # build blast database
    makeblastdb_cmd = f"makeblastdb -in {subject} -dbtype prot -out {subject}"
    p = Popen(makeblastdb_cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    blastp_cmd = f"blastp -query {query} -subject {subject} -num_threads {threads} -outfmt {outfmt} -out {out}"
    p = Popen(blastp_cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    return stdout, stderr

def high_identity_filter(blastp_out, identity=0.90, overlap=0.8):
    """
    Filter blastp result with identity > 0.90 and overlap > 80%
    """
    keep_set = set()
    remove_set = set()    
    gene_size = {}

    for line in blastp_out.split('\n'):
        line = line.strip()
        if line == "":
            continue
        items = line.split("\t")
        query_id = items[0]
        subject_id = items[1]
        subject_size = int(items[3])

        if query_id == subject_id:
            gene_size[query_id] = subject_size

    for line in blastp_out.split('\n'):
        line = line.strip()
        if line == "":
            continue
        items = line.split("\t")
        query_id = items[0]
        subject_id = items[1]
        subject_size = int(items[3])
        if query_id in remove_set:
            continue
        keep_set.add(query_id)
        if query_id == subject_id:
            continue
        
        identity = float(items[2])
        # identity > 90%
        if identity < 0.90:
            continue
        # overlap > 80%
        if subject_size < overlap * gene_size[query_id]:
            continue
        remove_set.add(subject_id)
    return keep_set,remove_set


def plot_exon_num_distribution(exon_num_list, upper_bound=19):
    from collections import Counter

    # Count occurrences of each exon number
    exon_num_counts = Counter(exon_num_list)

    # Modify counts to fit the specified categories
    exon_categories = {str(i): 0 for i in range(1, upper_bound+1)}
    exon_categories[f">{upper_bound}"] = 0
    for exon_num, count in exon_num_counts.items():
        if exon_num > upper_bound:
            exon_categories[f">{upper_bound}"] += count
        else:
            exon_categories[str(exon_num)] += count

    # Adjusting the bar chart to represent larger numbers in a compact form by capping the maximum length of the bar
    # Constants for ASCII chart
    MAX_BAR_LENGTH = 50  # Maximum length of the bar for the largest category
    largest_count = max(exon_categories.values())
    scale_factor = largest_count / MAX_BAR_LENGTH if largest_count > MAX_BAR_LENGTH else 1

    # Create adjusted ASCII bar chart
    adjusted_bar_chart = ""
    for category, count in exon_categories.items():
        adjusted_bar_length = int(count / scale_factor)  # Scale down the count to fit within the max bar length
        adjusted_bar_chart += f"{category}: " + "#" * adjusted_bar_length + "\n"

    print(adjusted_bar_chart)

def find_exon_num_threshold(exon_num_list):
    from collections import Counter

    # Count occurrences of each exon number
    exon_num_counts = Counter(exon_num_list)

    # Calculate the cumulative sum of gene counts
    cumulative_counts = []
    total_count = sum(exon_num_counts.values())
    cumulative_sum = 0
    for count in exon_num_counts.values():
        cumulative_sum += count
        cumulative_counts.append(cumulative_sum)

    # Find the threshold where the cumulative percentage is just over 80%
    threshold = None
    for i, cum_count in enumerate(cumulative_counts):
        if cum_count/total_count >= 0.8:
            threshold = i  # i corresponds to the index of sorted_exon_categories
            break

    # If a threshold is found, it corresponds to the category (exon count) that splits the data
    threshold_category = list(exon_num_counts.keys())[threshold] if threshold is not None else None

    threshold_category, threshold

def export(gene_dict, prefix):
    # step3: export the new gene_dict to gff3 and fasta
    print("export the new gene_dict to gff3 and fasta")
    gff3_out = f"{prefix}.gff3"

    with open(gff3_out, "w") as gff3_file:
        for gene in gene_dict.values():
            gff3_file.write(gene.to_gff3() + "\n")

    # step4: export the cds sequence to fasta
    cds_out = f"{prefix}.cds.fa"
    pep_out = f"{prefix}.pep.fa"
    
    with open(cds_out, "w") as cds_file:
        for gene in gene_dict.values():
            for tx in gene.transcripts:
                cds_file.write(f">{tx.id}\n{tx.orf.sequence}\n")

    with open(pep_out, "w") as pep_file:
        for gene in gene_dict.values():
            for tx in gene.transcripts:
                pep_file.write(f">{tx.id}\n{tx.orf.pep}\n")

def prepare_training_data(gff_file, fasta_file, gtf_source:str, prefix:str, min_exon_num: int, max_exon_num: int, min_orf_size ):
    """
    basic filtering: exon number < 2 or exon number > 1000
    """
    fasta = read_fasta(fasta_file)
    gene_dict = None

    if gtf_source == "stringtie":
        gene_dict = stringtie_loader(gff_file, fasta)

    # step1: parse gene_dict to filter transcripts with exon number < 2 or exon number > 1000
    gene_del_list = []
    for gene in gene_dict.values():
        if gene.transcripts is None:
            gene_del_list.append(gene.id)
            continue
        for tx in gene.transcripts:
            if tx.exon_num < min_exon_num or tx.exon_num > max_exon_num:
                gene.remove_transcript(tx.id)
        if gene.tx_num == 0:
            gene_del_list.append(gene.id)

    # delete genes with no transcripts
    for gene_id in gene_del_list:
        del gene_dict[gene_id]

    print(f"gene number: {len(gene_dict)}")

    # step2: select the gene with length > 3sigma - mean and < 3sigma + mean
    # step2.1: calculate the mean and sigma of gene length
    gene_len_list = []
    for gene in gene_dict.values():
        gene_len_list.append(gene.end - gene.start + 1)
    
    gene_len_list.sort()
    min_gene_size = gene_len_list[0]
    max_gene_size = gene_len_list[-1]
    mean, sigma = calculate_mean_and_sigma(gene_len_list)
    low_bound = mean - 3 * sigma
    high_bound = mean + 3 * sigma

    print(f"min_gene_size: {min_gene_size}")
    print(f"max_gene_size: {max_gene_size}")
    print(f"mean: {mean}")
    print(f"sigma: {sigma}")

    # step2.2: delete genes with length < low_bound or length > high_bound
    gene_del_list = []
    for gene in gene_dict.values():
        gene_size = gene.end - gene.start + 1
        if gene_size < low_bound or gene_size > high_bound:
            gene_del_list.append(gene.id)

    for gene_id in gene_del_list:
        del gene_dict[gene_id]

    print(f"gene number: {len(gene_dict)}")

    # step3: select the transcript with longest ORF for each gene and size > 300
    new_gene_dict :Dict[str, Gene] = {}
    for gene in gene_dict.values():
        seq = fasta[gene.chrom][gene.start-1:gene.end]
        if 'N' in seq or 'n' in seq: # skip the gene with N
            continue
        longest_orf = None
        longest_orf_len = 0
        longest_tx = None
        
        for tx in gene.transcripts:
            if tx.orf is not None and len(tx.orf.sequence) > longest_orf_len:
                longest_orf = tx.orf
                longest_orf_len = len(tx.orf.sequence)
                longest_tx = tx
            
        if longest_orf is not None and len(longest_orf.sequence) > min_orf_size :
            new_gene_dict[gene.id] = Gene(gene.id, gene.chrom, gene.start, gene.end, gene.strand)
            new_gene_dict[gene.id].add_transcript(longest_tx)

    # step4: self-alignment to remove redundant genes
    ## 输出临时的pep
    tmp_pep = f"tmp.pep.fa"
    with open(tmp_pep, "w") as tmp_pep_file:
        for gene in new_gene_dict.values():
            if gene.transcripts is None:
                continue
            for tx in gene.transcripts:
                tmp_pep_file.write(f">{tx.id}\n{tx.orf.pep}\n")

    ## 进行blastp
    stdout, stderr = diamond_blastp(tmp_pep, tmp_pep, outfmt=6)
    if len(stderr) > 0:
        print(stderr)
        exit(1)
    blastp_out = stdout.decode('utf-8')

    keep_set, _ = high_identity_filter(blastp_out)

    keep_set = {'.'.join(x.rsplit('.')[:-1]) for x in keep_set}

    # step5: remove genes with no transcripts and gene not in keep set
    gene_del_list = []
    for gene in new_gene_dict.values():
        if gene.transcripts is None:
            gene_del_list.append(gene.id)
            continue
        if gene.id not in keep_set:
            gene_del_list.append(gene.id)

    for gene in gene_del_list:
        del new_gene_dict[gene]

    # count the exon number of each gene
    exon_num_list = []
    for gene in new_gene_dict.values():
        for transcript in gene.transcripts:
            exon_num_list.append(transcript.exon_num)

    plot_exon_num_distribution(exon_num_list)
    # threshold_category, threshold = find_exon_num_threshold(exon_num_list)
    # print(f"threshold_category: {threshold_category}")
    # print(f"threshold: {threshold}")
    # output
    # save the new gene_dict to gff3 and fasta
    export(new_gene_dict, prefix)


def main(args):
    min_exon_num = args.min_exon_num
    max_exon_num = args.max_exon_num
    min_orf_size = args.min_orf_size
    prefix = args.prefix
    gff_file = args.gff
    fasta_file = args.ref
    gtf_source = args.format

    prepare_training_data(gff_file, fasta_file, gtf_source, prefix, min_exon_num, max_exon_num, min_orf_size )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='gene annotation adapter')
    parser.add_argument('-m', '--min_exon_num', type=int, default=2, help='min exon number')
    parser.add_argument('-M', '--max_exon_num', type=int, default=1000, help='max exon number')
    parser.add_argument('-s', '--min_orf_size', type=int, default=300, help='min orf size')
    parser.add_argument('-p','--prefix', type=str, help='prefix of output file', default='wgap')
    parser.add_argument('-f', '--format', type=str, help='format of input gff file', default='stringtie')
    parser.add_argument('ref', type=str, help='fasta file')
    parser.add_argument('gff', type=str, help='input gtf file')
    args = parser.parse_args()
    main(args)
