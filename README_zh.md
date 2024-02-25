# WGAP: Whole Genome Annotation Pipeline

Citation: [![DOI](https://zenodo.org/badge/363893963.svg)](https://zenodo.org/badge/latestdoi/363893963)

通常全基因组注释通常伴随着多个步骤，转录组数据组装，模型训练，从头预测等步骤。

BRAKER支持利用蛋白数据和转录组数据进行augustus的模型训练，最后输出augustus的预测结果，但是流程中使用到GeneMark在某些物种种表现不佳，导致augustus的模型训练存在问题。

MAKER是一个比较早的软件，在已有的基因从头预测模型，根据转录组和蛋白数据输出更加可靠的基因模型。但是该流程需要我们前期处理转录组数据，以及snap和augustus的模型需要自己一步步实现，比较繁琐。同时MAKER本质上还是从头预测软件，模型还是存在不少问题。

我从自己的实际需求出发，开发了一套流程，命名为WGAP (whole gneome annotation pipeline). 在这个流程中，我尽可能减少了对一些复杂软件的依赖，并且移除了最初对MAKER的依赖。不过还是需要一个EvidenceModeler, 只不过目前这个工具支持Singularity, 因此使用相对简单。

## 安装

WGAP的运行依赖于以下软件

- python=3.7
- pandas
- snakemake>=5.9.1
- click=7
- biopython=1.78
- augustus>=3.4.0
- stringtie>=2.1.7
- star
- minimap2
- miniprot
- samtools
- fastp
- gffread
- pblat>=2.5

手动安装

从GitHub上下载项目

```bash
git clone https://github.com/xuzhougeng/wgap
cd wgap
```

通过bioconda安装依赖环境

```bash
conda env create -n wgap -f=wgapenv.yml
```

安装wgap

```bash
pip install .
```

## 使用

假设你的输入如下

- 参考基因组: genome.fa
- RNA-seq组装得到GTF文件: tissue1.gtf, tissue2.gtf, ...


TODO: 增加RNA-seq的组装流程

### 使用EDTA预测TE

```Python
EDTA.pl -genome genome.fa -species others -step all -t 20
```

利用RepeatMasker对基因组进行mask

```Bash
RepeatMasker -xsmall -e rmblast -pa 20 -lib EDTA_TElib.fa  genome.fa -dir repeat
```

### 蛋白序列的回帖

为了速度，使用miniprot进行回帖 [TODO:加入自动化流程]

```Bash
mkdir homology && cd homology
miniprot -t 16 -d genome.mpi ../genome.fa
# 蛋白序列回帖
miniprot -Iut16 --gff genome.mpi homology.faa > homology.pep.gff
# （small exon 不一定找的对）
```

### 训练模型

二代测序数据质控，数据回帖，然后StringTie预测，得到初步的基因模型

```Python
# 准备模型的数据
wgap prepare gtf2train --gtf_source stringtie --prefix wgap tissue1.gtf genome.fa
```

训练模型

```Bash
mkdir -p model && cd model
wgap train augustus --prefix your_species --gene_number 1000 --intergenic_size 2000 wgap.gff3 ../genome.fa
# snap
wgap train snap --prefix snap --intergenic_size 1000 wgap.gff3 ../genome.fa
# glimmerhmm
wgap train glimmer --prefix glimmer --gene_number 2000 wgap.gff3 ../genome.fa

```

预测基因组：

```Python
# 先拆分
wgap predict split  --method te --repeat_masker ../repeat/genome.fa.out --gene_gff ../rna-seq/tissue1.gtf --homology ../homology/homology.pep.gff  ../genome.fa

# 在预测, 
wgap predict augustus --augustus_config_path $PWD/augustus/config/ --species your_species split_genome tmp
wgap predict snap --snap_hmm_file snap.hmm split_genome snap_predict
wgap predict glimmer --glimmerhmm_model_dir TrainGlimmM_glimmer split_genome glimmer_predict

```


### 整理EVM的输入

处理miniprot和rna-seq的输入文件

```Python
# rna-seq splice alignment
ls rna-seq/*gtf | while read gtf ; do wgap evm convert --source 1 $gtf  ;done

# protein splice alignment
wgap evm convert --source 2 -mm 2 homology/homology.pep.gff

# protein gene structure
wgap evm convert --source 2 -mm 1 homology/homology.pep.gff

```

获取权重文件

```Python
wgap evm init -ab "model/augustus.gff3 model/glimmer.gff3 model/snap.gff3" \
-tr "rna-seq/tissue1_evm.gff3 rna-seq/tissue2_evm.gff3 " \
-pr "homology/homology_spliced_alignment_evm.gff3" \
-ot homology/homology_gene_structure_evm.gff3
```

得到权重文件 weight.txt 和对应的输入

```Bash
ABINITIO_PREDICTION  augustus  1
ABINITIO_PREDICTION  glimmer  1
ABINITIO_PREDICTION  snap  1
TRANSCRIPT  tissue1  4
TRANSCRIPT  tissue2  4
PROTEIN  homology_gene_structure_evm  3
OTHER_PREDICTION  homology_spliced_alignment_evm  5
```

使用EVM

```R
singularity exec -e  ~/images/EVidenceModeler.v2.1.0.simg EVidenceModeler \
    --sample_id your_species \
    --genome genome.fa \
    --weights weights.txt \
    --CPU 64 \
    --gene_predictions gene_predictions.gff3 \
    --protein_alignments protein_alignments.gff3 \
    --transcript_alignments transcript_alignments.gff3 \
    --segmentSize 100000 \
    --overlapSize 10000
```

输出结果，可以用ingenannot或者PASApipeline增加UTR