# WGAP: Whole Genome Annotation Pipeline

Citation: [![DOI](https://zenodo.org/badge/363893963.svg)](https://zenodo.org/badge/latestdoi/363893963)

通常全基因组注释通常伴随着多个步骤，转录组数据组装，模型训练，从头预测等步骤。

BRAKER支持利用蛋白数据和转录组数据进行augustus的模型训练，最后输出augustus的预测结果，但是流程中使用到GeneMark在某些物种种表现不佳，导致augustus的模型训练存在问题。

MAKER是一个比较早的软件，在已有的基因从头预测模型，根据转录组和蛋白数据输出更加可靠的基因模型。但是该流程需要我们前期处理转录组数据，以及snap和augustus的模型需要自己一步步实现，比较繁琐。

为了简化基因组的注释流程，我基于snakemake开发了WGAP流程。WGAP流程的主体是MAKER，通过maker基于转录组和蛋白数据输出较为可靠的模型，之后基于模型对从头预测软件进行迭代训练，最后输出基因注释结果。

WGAP能够同时使用二代和三代转录组数据，提高注释结果（三代转录组处理流程还在优化中）,支持从Swiss/Prot中下载相关物种的高质量蛋白序列。

WGAP能够自动处理多种场景

- 仅有转录组数据，进行模型训练和基因预测
- 仅有蛋白数据，进行模型训练和基因预测
- 已有模型，基于转录组和蛋白序列进行模型优化然后进行基因预测
- 已有模型和蛋白数据，不需要预测，直接输出结果


## 安装

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
python setup.py install
```

直接使用bioconda安装(TODO, 等待软件稳定)

## 使用方法

建立一个项目文件夹，例如annotation

```bash
mkdir -p annotation && cd annotation 
```

初始化配置文件， 这一步会复制config.yaml 和 samples.csv 文件到annotation目录下

```bash
wgap init
```

下载物种的蛋白序列, 以植物为例

```bash
wgap download -s plants -d sprot
```

注, 这一步可能会因为网络原因下载失败，多试几次即可。

接着修改config.yaml中的配置信息

必须填写的内容

- genome: 基因组序列的文件路径
- specie_name: 物种名，后续的augustus的模型
- protein: 蛋白序列的位置

对于转录组数据, 本流程使用fastp, star, stringtie流程获取组装的转录本序列

- transcript_assemble: True 或者 False, 是否组装转录本, 如果不需要组装，那么sample一行会被忽略
- sample: "samples.csv"
- transcript_gtf: 已有的GFF位置信息
- transcript_fasta: 已有的组装的转录本数据

对于 transcript_gtf 文件，maker要求gtf文件的第二列是match和math_part形式，可以通过如下命令进行转换

```bash
gtf=你的输入gtf文件
gffread -E $gtf -o- | sed "s#transcript#match#g" |sed "s#exon#match_part#g" > rna_seq.gff3
```

其中, samples.csv一共分为5列，其中第二列指的是tissue, 相同tissue的多组重复会被当作一个。
（目前暂时只支持二代转录组数据，后续会增加三代转录组的支持）

基因训练相关参数

- training_model: True 或者 False, 当为False时，则snap和augustus模型不会训练，同时必须保证snaphmm和augustus_species其中一个不为空
- training_snap_round: 12, 表示训练2轮，如果填写123，则表示训练三轮
- training_augustus_round: 2, 表示在第2轮种训练augustus, 如果填写123, 表示和snap一样训练三轮

已有的snap和augustus模型

- snaphmm:  hmm路径
- augustus_species: 物种名

当已有snap和augustus模型时，wgap会使用已有的模型的进行优化。

运行命令

```bash
export LIBDIR=~/miniconda3/envs/wgap/share/RepeatMasker/Libraries
wgap run all -j 80 
```

使用bioconda安装的RepeatMasker需要配置LIBDIR环境变量，否则会报错。

## 输出结果

运行结束之后，关键的输出文件是如下三个

- maker.gff: maker的注释结果
- abinit.gff: snap, augustus的从头注释结果中和maker不重叠的基因
- other.gff: 记录蛋白,转录组,重复序列位置信息

maker的输出结果相对于其他注释软件会比较少，意味着有些从头预测的基因由于缺少证据被maker所过滤掉。

为了提高基因数目，我开发 `rescue` 模块，用来从 abinit.gff中恢复一些较为可信的基因。

我们可以使用InterProScan对maker目录下的maker_annotation.all.maker.non_overlapping_ab_initio.proteins.fasta
进行功能域预测，保留具有功能域的基因。

```bash
interproscan.sh -cpu 80 -appl pfam -dp -f TSV  -t p -i maker/maker.all.maker.non_overlapping_ab_initio.proteins.fasta -o output.iprscan &
cut -f 1 output.iprscan | sort -u > resue_list.txt
```

根据基因列表，将后续基因增加到maker.gff中, 默认输出结果是 maker_recue.gff

```bash
wgap rescue maker.gff abinit.gff resue_list.tx
```

maker的默认命名格式不太适合用于后续展示，我开发了 `rename` 模块，用于对结果进行重命名

```bash
wgap rename -p AT maker_recue.gff annotation.gff
```

输出结果形如 AT1G00010, AT1G00020 其中AT是物种拉丁名缩写, 1G表示1号染色体, 00010, 00020分别表示第一个基因, 第二个基因。
对于非染色体上的基因，例如,则直接以AT000010的形式进行编号

> 目前不支持线粒体和叶绿体的编号，后续会尽量支持更多类型的命名形式

## 手工校正

由于目前的注释流程都无法保证百分百的准确率，意味着后续还需要进行手工调整基因结构。为了方面对原来的gff文件进行调整，我开发 `update` 模块。

```bash
wgap -o maker_update.gff maker_recue.gff apollo_exported.gff
```

不过目前该模块要求的输入 gff 文件的第9列必须包括 description 信息，根据 description 中的内容对原来基因模型进行调整

- 基因名: 
  - 基因名在原本物种中不存在时，新增基因，
  - 基因名在原本物种中存在时, 修改基因模型
- 基因名delete: 删除该基因


调整完之后，还可以继续用 `rename` 进行重命名。

## 已知问题:

运行到一半的时候，会提示任务被挂起

```bash
suspended (tty input)  wgap run -j 80
```

解决方法，使用 `bg %wgap` 即可让任务在后台继续运行。