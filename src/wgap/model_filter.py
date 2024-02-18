import pandas as pd
import pyranges as pr
import argparse

def filter_gff_by_rank(model_rank_csv, gff_file, output_gff_file):
    
    # 读取 model rank 数据
    model_rank_df = pd.read_csv(model_rank_csv, sep="\t")

    gene_pr = pr.read_gff3(gff_file)

    # 获取 rank 大于 1 的行的索引
    gene_id = model_rank_df[model_rank_df['rank'] > 1]['index'].to_list() 
    
    # 找到mRNA的行
    mrna_pr = gene_pr[gene_pr.Feature == 'mRNA']
    mrna_id = mrna_pr.ID[mrna_pr.Parent.isin(gene_id)].to_list()

    # 过滤:
    filtered_pr = gene_pr[gene_pr.ID.isin(gene_id) |  gene_pr.Parent.isin(gene_id) | gene_pr.Parent.isin(mrna_id)]

    filtered_pr.to_gff3(output_gff_file)
    

def main():
    parser = argparse.ArgumentParser(description="过滤 GFF 文件基于 model rank.")
    parser.add_argument("model_rank_csv", help="model rank CSV 文件路径")
    parser.add_argument("gff_file", help="输入 GFF 文件路径")
    parser.add_argument("output_gff_file", help="输出 GFF 文件路径")

    args = parser.parse_args()

    filter_gff_by_rank(args.model_rank_csv, args.gff_file, args.output_gff_file)
    print(f"Filtered GFF data saved to {args.output_gff_file}")

if __name__ == "__main__":
    main()
