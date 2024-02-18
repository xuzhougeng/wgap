import sys

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


def main(gff_file,output_file):
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


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <GFF3 file> <ouput file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
