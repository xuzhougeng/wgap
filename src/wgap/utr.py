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


def update_gene(gene: Gene, tx_list: List[Transcript]) -> Gene:
    """
    Process a gene to find the transcript with the maximum overlap ratio.
    
    :param gene: A Gene object that includes chrom, start, end, and strand information.
    :param tx_list: A list of Transcript objects.
    :return: A Gene object with the update transcript 
    """

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


def add_utr_to_gene_model(gff_file, gtf_file_list, out_gff_file, thread_num = 20):

    # Load EVM gff3 file into
    gene_dict:Dict[str, Gene] = gff3_loader(gff_file)
    
    # query ranges of all genes
    query_ranges = pr.from_dict({
        'Chromosome': [gene.chrom for gene in gene_dict.values()],
        'Start': [gene.start for gene in gene_dict.values()],
        'End': [gene.end for gene in gene_dict.values()],
        'Strand': [gene.strand for gene in gene_dict.values()],
        'GeneID': list(gene_dict.keys())
    })

    # Load gtf files into PyRanges
    gtf_list = []
    for gtf_file in gtf_file_list:
        gtf_source = basename(gtf_file).split('.')[0]
        gtf_pr = pr.read_gtf(gtf_file)
        gtf_pr.Source = gtf_source
        gtf_list.append(gtf_pr)

    # merge gtf into one PyRanges
    gtf_pr = pr.concat(gtf_list)

    # 筛选出转录本和外显子
    tx_gtf_pr = gtf_pr[gtf_pr.Feature == 'transcript']
    tx_gtf_pr = tx_gtf_pr.drop(['Score', 'Frame', 'exon_number'])
    exon_gtf_pr = gtf_pr[gtf_pr.Feature == 'exon']

    exon_gtf_pr.exon_id = exon_gtf_pr.transcript_id + '.exon.' + exon_gtf_pr.exon_number.astype(str)
    exon_gtf_pr = exon_gtf_pr.drop([ 'Feature', 'Source', 'Score', 'Strand', 'Frame', 'gene_id', 'FPKM', 'TPM', 'cov','exon_number'])

    # 合并query_ranges 和 tx_gtf_pr, 得到每个gene的转录本信息
    query_tx_pr = query_ranges.join(tx_gtf_pr, suffix= '_tx')

    full_df = query_tx_pr.df.merge(exon_gtf_pr.df, on='transcript_id', suffixes=('', '_exon'))

    # 首先按 GeneID 和 transcript_id 分组
    grouped = full_df.groupby(['GeneID', 'transcript_id'])

    def create_transcript(group):
        """根据外显子信息创建转录本对象"""
        exon_list = [Exon(row['exon_id'], row['Chromosome_exon'], row['Start_exon'], row['End_exon'], row['Strand']) for idx, row in group.iterrows()]
        first_row = group.iloc[0]
        return Transcript(first_row['transcript_id'], first_row['Chromosome'], first_row['Start_tx'], first_row['End_tx'], first_row['Strand_tx'], exon_list)

    # 使用 apply() 通过每个组（每个基因的每个转录本）创建转录本对象
    transcripts = grouped.apply(create_transcript)

    # 重新按 GeneID 分组以聚合每个基因的转录本列表
    gene_to_tx_list = transcripts.groupby(level=0).apply(list).to_dict()

    # update gene model for each gene in gene_dict
    from concurrent.futures import ProcessPoolExecutor

    def parallel_update_gene(gene_id: str, gene: Gene, tx_list: List[Transcript]) -> Gene:
        """封装 update_gene 函数以适应并行调用的需求"""
        if len(tx_list) == 0:
            return gene_id, gene
        updated_gene = update_gene(gene, tx_list)
        return gene_id, updated_gene

    # 存储更新后的基因模型
    updated_gene_dict = {}

    with ProcessPoolExecutor(max_workers=thread_num) as executor:
        futures = []
        for gene_id, gene in gene_dict.items():
            # 提交并行任务
            tx_list = gene_to_tx_list.get(gene_id, [])
            future = executor.submit(parallel_update_gene, gene_id, gene, tx_list)
            futures.append(future)
        
        # 收集并行执行的结果
        for future in futures:
            gene_id, updated_gene = future.result()
            updated_gene_dict[gene_id] = updated_gene

    return updated_gene_dict

    # # save to gff3
    # with open(out_gff_file, 'w') as f:
    #     for gene in updated_gene_dict.values():
    #         f.write(gene.to_gff3())
    #         f.write('\n')