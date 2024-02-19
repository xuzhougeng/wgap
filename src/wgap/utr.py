from typing import List,Dict, Tuple, Optional
import concurrent.futures
from functools import partial
from os.path import basename

import pyranges as pr
from pyranges import PyRanges

from .model_loader import gff3_loader
from .gene import Exon, Transcript, Gene


def calculate_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Calculate the overlap length between two intervals.

    :param a_start: The start position of interval A
    :param a_end: The end position of interval A
    :param b_start: The start position of interval B
    :param b_end: The end position of interval B
    :return: The overlap length between the two intervals
    """

    return max(0, min(a_end, b_end) - max(a_start, b_start))

def structural_similarity(transcript_a: List[Transcript], transcript_b: List[Transcript]) -> int:

    """Calculate the structural similarity (or structural distance) between two transcripts.
    :param transcript_a: The first transcript
    :param transcript_b: The second transcript
    :return: The structural distance between the two transcripts
    """

    overlap_total = 0
    for exon_a in transcript_a.exons:
        for exon_b in transcript_b.exons:
            overlap_total += calculate_overlap(exon_a.start, exon_a.end, exon_b.start, exon_b.end)
    
    # suppose the length of each exon is calculated and stored in the transcript object
    total_length_a = sum([exon.end - exon.start for exon in transcript_a.exons])
    total_length_b = sum([exon.end - exon.start for exon in transcript_b.exons])
    
    # structural distance
    structural_distance = (total_length_a + total_length_b) - 2 * overlap_total
    return structural_distance

def find_most_similar_transcript(transcript_a : Transcript, b_list: List[Transcript]) -> Optional[Transcript]:
    """Find the transcript in B_list that is most structurally similar to A.

    :param transcript_a: The first transcript
    :param b_list: List of candidate transcripts
    :return: The transcript most structurally similar to A

    """

    most_similar = None
    smallest_distance = float('inf')
    for transcript_b in b_list:
        current_distance = structural_similarity(transcript_a, transcript_b)
        if current_distance < smallest_distance:
            smallest_distance = current_distance
            most_similar = transcript_b
    return most_similar

def identify_utr(exon_list_A: List[Exon], exon_list_B: List[Exon], strand: str) -> Tuple[List[Exon], List[Exon]]:
    """
    Identify the 5'UTR and 3'UTR regions based on two sets of exons and the strand.

    :param exon_list_A: The first set of exons represented as Exon objects.
    :param exon_list_B: The second set of exons represented as Exon objects.
    :param strand: The strand of the gene ("+" or "-").
    :return: A tuple containing two lists of Exon objects representing the 5'UTR and 3'UTR regions.
    """

    if strand == '+':
        lowest_start_A = min(exon.start for exon in exon_list_A)
        highest_end_A = max(exon.end for exon in exon_list_A)
        utr_5 = [exon for exon in exon_list_B if exon.start < lowest_start_A]
        utr_3 = [exon for exon in exon_list_B if exon.end > highest_end_A]
    elif strand == '-':  
        highest_start_A = max(exon.start for exon in exon_list_A)
        lowest_end_A = min(exon.end for exon in exon_list_A)
        utr_5 = [exon for exon in exon_list_B if exon.end > highest_start_A]
        utr_3 = [exon for exon in exon_list_B if exon.start < lowest_end_A]
    else:
        raise ValueError(f"Invalid strand: {strand}, only supprt '+' or '-'.")
    
    return utr_5, utr_3

def merge_exons_and_utrs(exon_list: List[Exon], utr_5:List[Exon], utr_3:List[Exon]):
    """
    
    :param exon_list: A list of Exon objects.
    :param utr_5: A list of Exon objects representing the 5'UTR region.
    :param utr_3: A list of Exon objects representing the 3'UTR region.
    :return: A list of Exon objects representing the merged exons and UTR regions.
    
    """
    # 假设 exon_list, utr_5, utr_3 都是 Exon 对象的列表
    all_exons = sorted(exon_list + utr_5 + utr_3, key=lambda x: x.start)

    merged_exons = []
    for current in all_exons:
        if not merged_exons:
            merged_exons.append(current)
        else:
            last = merged_exons[-1]
            if current.start <= last.end:  # 如果当前exon与上一个exon重叠
                # 更新重叠exon的结束位置为当前和上一个exon结束位置的最大值
                last.end = max(last.end, current.end)
                # 注意：这里假设Exon对象允许直接修改end属性
            else:
                merged_exons.append(current)

    return merged_exons

def update_gene(gene: Gene, gtf_pr: PyRanges) -> Gene:
    """
    Process a gene to find the transcript with the maximum overlap ratio.
    
    :param gene: A Gene object that includes chrom, start, end, and strand information.
    :param gtf_pr: A PyRanges object containing gene annotations .
    :return: A Gene object with the update transcript 
    """

    query_range = pr.PyRanges(chromosomes=[gene.chrom], starts=[gene.start], ends=[gene.end], strands=[gene.strand])

    tx_gtf_pr = tx_gtf_pr = gtf_pr[gtf_pr.Feature == 'transcript'] # gtf_pr with only transcript
    overlapping_regions = tx_gtf_pr.overlap(query_range)

    if len(overlapping_regions) == 0:
        return gene

    tx_list = []
    for _,row in overlapping_regions.df.iterrows():
        # 直接使用DataFrame的能力来过滤和处理
        exons_df = gtf_pr[ (gtf_pr.Feature== 'exon') & (gtf_pr.transcript_id== row['transcript_id']) ].df
        exons_df['exon_id'] = exons_df['transcript_id'] + '.exon.' + exons_df['exon_number'].astype(str)
        exon_list = [Exon(row['exon_id'], row['Chromosome'], row['Start'], row['End'], row['Strand']) for _, row in exons_df.iterrows()]

        tx = Transcript(row['transcript_id'], row['Chromosome'], row['Start'], row['End'], row['Strand'], exon_list)
        tx_list.append(tx)


    qry_tx = gene.transcripts[0]

    # 寻找最大重叠比例的Transcript
    best_tx = find_most_similar_transcript(qry_tx, tx_list)

    #  如果best_tx的exon总长度是qry_tx的两倍以上，则输出警告信息
    if best_tx is not None:
        best_tx_exon_length = sum(e.end - e.start + 1 for e in best_tx.exons)
        qry_tx_exon_length = sum(e.end - e.start + 1 for e in qry_tx.exons)
        if best_tx_exon_length > 2 * qry_tx_exon_length:
            # logging 
            print(f"Warning: The exon length of the best matching transcript {best_tx.id} is more than twice the length of the query transcript {qry_tx.id}.")

            return gene
    
    # 识别UTR
    utr_5, utr_3 = identify_utr(qry_tx.exons, best_tx.exons, gene.strand)
    
    # 合并exons和UTR
    updated_exons = merge_exons_and_utrs(qry_tx.exons, utr_5, utr_3)

    # 新建Gene和Transcript对象
    update_start = min(exon.start for exon in updated_exons)
    update_end = max(exon.end for exon in updated_exons)

    update_tx = Transcript(qry_tx.id, qry_tx.chrom, update_start, update_end, qry_tx.strand, updated_exons)
    # add cds to transcript
    for cds in qry_tx.cds:
        update_tx.add_cds(cds)
    update_gene = Gene(gene.id, gene.chrom, update_start, update_end, gene.strand, [update_tx])

    return update_gene


def update_gene_wrapper(gene_item, gtf_pr):
    """Simple wrapper function to update_gene with a single argument.
    
    """
    gene_id, gene = gene_item
    return gene_id, update_gene(gene, gtf_pr)


def add_utr_to_gene_model(gff_file, gtf_file_list, thread_num = 10):
    """
    
    """
    
    # Load EVM gff3 file into
    gene_dict:Dict[str, Gene] = gff3_loader(gff_file)

    # Load gtf files into PyRanges
    gtf_list = []
    for gtf_file in gtf_file_list:
        gtf_source = basename(gtf_file).split('.')[0]
        gtf_pr = pr.read_gtf(gtf_file)
        gtf_pr.Source = gtf_source
        gtf_list.append(gtf_pr)

    # merge gtf
    gtf_pr = pr.concat(gtf_list)

    process_with_gtf = partial(update_gene_wrapper, gtf_pr=gtf_pr)

    # 使用ThreadPoolExecutor并发处理每个基因
    processed_genes = {}

    #  with concurrent.futures.ThreadPoolExecutor() as executor:
    with concurrent.futures.ProcessPoolExecutor( max_workers= thread_num ) as executor:
        # 使用map函数处理所有基因
        results = executor.map(process_with_gtf, gene_dict.items())

        # 从结果中收集处理后的基因数据
        for gene_id, processed_gene in results:
            processed_genes[gene_id] = processed_gene

    
    return processed_genes

