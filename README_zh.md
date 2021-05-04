# WGAP: Whole Genome Annotation Pipeline

WGAP是一个基于MAKER的全基因组基因结构注释流程，

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

## 使用

建立一个项目文件夹，例如annotation

```bash
mkdir -p annotation && cd annotation 
```

初始化配置文件

```bash
wgap init

```

下载物种的蛋白序列, 以植物为例

```bash
wgap download -s plants -d sprot
```

可能会因为网络原因下载失败，多试几次即可。

接着修改config.yaml中的配置信息

必须填写的内容

- genome: 基因组序列的文件路径
- specie_name: 物种名，后续的augustus的模型
- protein: 蛋白序列的位置

对于转录组数据, 本流程使用fastp, star, stringtie流程获取组装的转录本序列

- transcript_assemble: True 或者 False, 是否组装转录本, 如果不需要组装，那么sample一行会被忽略
- sample: "samples.csv"
- gtf: 已有的转录本数据

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

引用方式

[![DOI](https://zenodo.org/badge/363893963.svg)](https://zenodo.org/badge/latestdoi/363893963)