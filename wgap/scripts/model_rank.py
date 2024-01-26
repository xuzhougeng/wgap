import pyarrow as pr
import pandas as pd
from typing import List
from joblib import Parallel, delayed
import multiprocessing


def read_STAR_SJ_out_tab(file: str) -> pd.DataFrame:
    """Read STAR SJ.out.tab file
    
    Args:
        file (str): STAR SJ.out.tab file
    
    Returns:
        pr.Table: STAR SJ.out.tab file
    """

    # Moitf:  Intron types, where 0 indicates none of the following functions, 1 represents GT/AG, 2 represents GT/AC, 3 represents GC/AG, 4 represents GT/GC, 5 represents AT/GC, 6 represents GT/AT
    df = pd.read_csv(file, sep="\t", names=["Chromosome", "Start", "End", "Strand", "Motif", "Annotated", "Unique", "Multi", "MaxOverhang"])
    df.Start = df.Start - 1
    df.Strand = df.Strand.map({1:"+", 2:"-"})

    return df

def create_hints_from_STAR_SJ_out_tab(file_list: List[str]) -> List[str]:
    """Create hints from STAR SJ.out.tab file
    
    Args:
        file_list (List[str]): STAR SJ.out.tab file list
    
    Returns:
        pd.DataFrame: hints
    """
    # read STAR SJ.out.tab file and concat them

    df = pd.concat([read_STAR_SJ_out_tab(file) for file in file_list])

    # add SJ_id column: 
    df["SJ_id"] = df['Chromosome'].astype(str) + ":" + df['Start'].astype(str) + "-" + df['End'].astype(str) 

    # sort by Chromosome and Start
    df = df.sort_values(by=["SJ_id"])

    # 定义一个函数来处理每个分组
    def filter_group(group):
        # 如果 Strand 和 Motif 在组内不是唯一的，则返回 MaxOverHang 最大的记录
        if group['Strand'].nunique() == 1 and group['Motif'].nunique() == 1:
            sorted_group = group.sort_values(by=['MaxOverhang', 'Unique', 'Multi'], ascending=False)
            return sorted_group.iloc[0]
        # 否则，返回整个组
        return group

    def applyParallel(dfGrouped, func):

        max_cores = multiprocessing.cpu_count()
        # use 1/4 cores to run at lease 4
        n_jobs = max(max_cores//4, 4)

        retLst = Parallel(n_jobs= n_jobs )(delayed(func)(group) for name, group in dfGrouped)
        return retLst

    res = applyParallel(df.groupby(['SJ_id']), filter_group)

    df = pd.DataFrame(res)

    # 过滤掉 Uniq + Multi < 10 
    df = df[df.Unique + df.Multi >= 10]

    hints = df.Chromosome.astype(str) + ":" + df.Start.astype(str) + "-" + df.End.astype(str) + ":" + df.Strand.astype(str)

    return hints.to_list()

def homology_match(df):
    # 用于评估基因坐标与同源蛋白的重叠率
    # =（完全相等), <(被包含) , >(包含), x（没有证据), p5(overlap > 0.5), mX（表示存在多个基因）
    if df.shape[0] == 1:
        if ((df['Start'] == df['Start_b']) & (df['End'] == df['End_b'])).all():
            return "="
        elif ((df['Start'] > df['Start_b']) & (df['End'] < df['End_b'])).all():
            return "<"
        elif ((df['Start'] < df['Start_b']) & (df['End'] > df['End_b'])).all():
            return ">"
        else:
            # calc the overlap rati
            overlap = min(df['End'].iloc[0], df['End_b'].iloc[0]) - max(df['Start'].iloc[0], df['Start_b'].iloc[0])
            overlap_ratio = overlap / (df['End'].iloc[0] - df['Start'].iloc[0])
            if overlap_ratio >= 0.5:
                return "p"
            else:
                return "x"
    else:
        return f"m{df.shape[0]}"

def model_rank(intron_support_ratio:float, homology_overlap:str):
    intron_rank = 0
    if intron_support_ratio == 1:
        intron_rank  = 2
    elif intron_support_ratio == 0:
        intron_rank = 0
    else:
        intron_rank = 1

    homology_rank = 0
    if homology_overlap == "x":
        homology_rank = 0
    elif homology_overlap.startswith("m"):
        homology_rank = 1
    elif homology_overlap == "=":
        homology_rank = 3
    else:
        homology_rank = 2

    return intron_rank + homology_rank

def main(augustus_gff, homology_gff: str, SJ_out_list: List[str] ):

    #  基本思路
    ## 1. 从augutus的gff文件得到每个基因位置
    ## 2. 从STAR的SJ.out.tab文件得到所有的可变剪接位点
    ## 3. 从miniprot的gff文件得到所有同源蛋白的位置
    ## 4. 基于STAR和miniprot对结果进行打分
    
    gff = pr.read_gff3(augustus_gff)
    
    gene_list = gff[gff.Feature == "gene"].ID.to_list()
    gene_dict = {}
    for gene in gene_list:
        gene_dict[gene] = { "intron_support_ratio": 0, "homology_overlap": "x" }

    # 从STAR的SJ.out.tab文件得到所有的可变剪接位点
    hints = create_hints_from_STAR_SJ_out_tab(SJ_out_list)

    # 获取内含子位置
    ## 提取 CDS 区域
    cds = gff[gff.Feature == "CDS"]
    ## 提取基因区域
    genes = gff[gff.Feature == "gene"]
    ## 反选内含子
    introns = genes.subtract(cds)
    
    ## 计算内含子支持率: 内含子存在于hints中的数量 / 内含子总数
    introns.Coord = introns.Chromosome.astype(str) + ":" + introns.Start.astype(str) + "-" + introns.End.astype(str) + ":" + introns.Strand.astype(str)
    introns.In_Hints = introns.Coord.isin(hints)
    support_ratio = introns.df.groupby(["ID"]).In_Hints.sum() / introns.df.groupby(["ID"]).In_Hints.count()

    for index, value in support_ratio.items():
        # 'index' is the index of the item in the series
        # 'value' is the value of the item in the series
        gene_dict[index]["intron_support_ratio"] = value

    # 同源蛋白
    homology_pr = pr.read_gff3( homology_gff )
    homology_pr = homology_pr[homology_pr.Rank == '1' ]
    
    # 基于rank筛选
    overlaps_pr = gff[gff.Feature == "mRNA"].join(homology_pr[homology_pr.Feature == "mRNA"])
    overlaps_df = overlaps_pr.df

    homo_evidence = overlaps_df.groupby(['Parent']).apply(lambda x: homology_match(x))
    for index, value in homo_evidence.items():
        # 'index' is the index of the item in the series
        # 'value' is the value of the item in the series
        gene_dict[index]["homology_overlap"] = value

    # rank the gene_dict
    for gene in gene_dict.keys():
        gene_dict[gene] = model_rank(**gene_dict[gene])

    df = pd.DataFrame.from_dict(gene_dict, orient='index').reset_index()
    df['model_rank'] = df.apply(lambda x: model_rank(x['intron_support_ratio'], x['homology_overlap']), axis=1)

    # 保存结果
    df.round(2).to_csv("gene_model_rank.csv", sep="\t", index=None)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="gene model rank")

    # 添加参数
    parser.add_argument("augustus_gff", help="路径到 Augustus GFF 文件")
    parser.add_argument("homology_gff", help="路径到同源性 GFF 文件")
    parser.add_argument("SJ_out_list", nargs='+', help="一系列 SJ 输出文件的路径")

    # 解析参数
    args = parser.parse_args()

    # 调用 main 函数
    main(args.augustus_gff, args.homology_gff, args.SJ_out_list)