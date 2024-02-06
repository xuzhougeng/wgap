from wgap.scripts.seq_utils import read_fasta, reverse_complement
from wgap.scripts.gene import gff3_loader
import random
import subprocess
import gzip
import sys
import os
import shutil
import re

def run_command(cmd, stdout_file=None, stderr_file=None):
    # 使用subprocess.run执行命令
    try:
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # 打印输出结果
        if stdout_file is not None:
            with open(stdout_file, "w") as fh:
                fh.write(result.stdout.decode())
        if stderr_file is not None:
            with open(stderr_file, "w") as fh:
                fh.write(result.stderr.decode())
    except subprocess.CalledProcessError as e:
        print("Error executing command:", e)

def mirror_intervals_within_range(A, sub_intervals):
    """
    将子区间列表中的每个区间翻转到主区间中。
    """
    # 计算主区间的长度
    # A_length = A[1] - A[0]
    
    # 初始化结果列表
    reversed_intervals = []
    
    # 遍历子区间列表
    for interval in sub_intervals:
        # 计算每个子区间的开始和结束点相对于主区间的新位置
        start_relative = interval[0] - A[0]
        end_relative = interval[1] - A[0]
        
        # 计算翻转后的位置
        new_start = A[1] - end_relative
        new_end = A[1] - start_relative
        
        # 将翻转后的子区间添加到结果列表中
        reversed_intervals.append((new_start, new_end))
    
    # 返回翻转后的子区间列表
    return reversed_intervals

def prepare_snap_training_data(fasta_path, gff_path, zff_path):
    """
    准备并打印基于FASTA和GFF文件的SNAP训练数据。

    参数:
    fasta_path -- FASTA文件的路径
    gff_path -- GFF文件的路径
    """
    gene_exons = {}
    try:
        gff_file = gzip.open(gff_path, 'rt') if gff_path.endswith('.gz') else open(gff_path, 'r')
        with gff_file as file:
            for line in file:
                if line.startswith('#'):
                    continue
                chromosome, _, feature, start, end, _, strand, _, attributes = line.strip().split('\t')
                if feature != 'CDS':
                    continue
                gene_name = [attr for attr in attributes.split(';') if attr.startswith('Parent=')][0].split('=')[1]
                start, end = (end, start) if strand == '-' else (start, end)
                if chromosome not in gene_exons:
                    gene_exons[chromosome] = {}
                if gene_name not in gene_exons[chromosome]:
                    gene_exons[chromosome][gene_name] = []
                gene_exons[chromosome][gene_name].append((start, end))
    except FileNotFoundError:
        sys.exit(f"File not found: {gff_path}")

    sequence_ids = []
    try:
        fasta_file = gzip.open(fasta_path, 'rt') if fasta_path.endswith('.gz') else open(fasta_path, 'r')
        with fasta_file as file:
            for line in file:
                if line.startswith('>'):
                    sequence_id = line.strip().split('>')[1]
                    sequence_ids.append(sequence_id)
    except FileNotFoundError:
        sys.exit(f"File not found: {fasta_path}")

    # 打印结果
    with open(zff_path, 'w') as file:
        for chromosome in sequence_ids:
            file.write(f">{chromosome}\n")
            if chromosome in gene_exons:
                for gene in gene_exons[chromosome]:
                    for exon in gene_exons[chromosome][gene]:
                        file.write(f"Exon {exon[0]} {exon[1]} {gene}\n")

def train_snap_model(gff3_file, genome_fasta, prefix = "snap", intergenic_size=1000):

    # create a directory to store the training data
    train_dir = "snap_train"
    if not os.path.exists(train_dir):
        os.makedirs(train_dir)
    
    # get real path of gff3_file and genome_fasta
    gff3_file = os.path.realpath(gff3_file)
    genome_fasta = os.path.realpath(genome_fasta)
    
    zff_file = os.path.realpath(f"{prefix}.zff")
    prepare_snap_training_data(genome_fasta, gff3_file, zff_file)
    
    # split genes into various categories
    cmd = f"cd {train_dir} && fathom -categorize {intergenic_size} {zff_file} {genome_fasta}"
    run_command(cmd)

    #  export all of the uni genes into their plus-stranded versions
    cmd = f"cd {train_dir} && fathom -export {intergenic_size} -plus uni.ann uni.dna"
    run_command(cmd)

    # valite again
    # cmd = f"cd {train_dir} && fathom -validate export.ann export.dna"
    # run_command(cmd)

    # create a large number of model files.
    cmd = f"cd {train_dir} && forge export.ann export.dna"
    run_command(cmd)

    cmd = f"hmm-assembler.pl {prefix} {train_dir}  > {prefix}.hmm"
    run_command(cmd)

def prepare_glimmer_training_data(gff3_file, genome_fasta, gene_number=500, prefix="glimmer"):
    fasta_dict = read_fasta(genome_fasta)

    gene_models = gff3_loader(gff3_file, fasta_dict)
    # random select 1000 genes from gene_models
    if len(gene_models) >= gene_number:
        # 将字典的项转换为列表，然后从中随机选择100个项
        random_items = random.sample(list(gene_models.items()), gene_number)
        # 转换回字典
        gene_models = dict(random_items)

    mfasta_file = f"{prefix}.mfasta"
    exon_file = f"{prefix}.exon_coords"

    with open(mfasta_file, "w") as fh1, open(exon_file, "w") as fh2:
        for _, gene in gene_models.items():
            for tx in gene.transcripts:
                full_seq = fasta_dict[tx.chrom][tx.start:tx.end]
                if tx.strand == "-":
                    full_seq = reverse_complement(full_seq)
                fh1.write(f">{tx.id}\n{full_seq}\n")
                
                cds_list = []
                for cds in tx.cds:
                    cds_list.append([cds.start, cds.end])

                # reorder cds_list

                cds_list.sort(key=lambda x: x[0])
                if tx.strand == "-":
                    cds_list = mirror_intervals_within_range([tx.start, tx.end], cds_list)
                    cds_list.sort(key=lambda x: x[0])

                for cds in cds_list:
                    cds_start = cds[0] - tx.start + 1
                    cds_end = cds[1] - tx.start
                    fh2.write(f"{tx.id}\t{cds_start}\t{cds_end}\n")
                    # print(cds_start ,cds_end)
            fh2.write(f"\n")
                
def train_glimmer_model(mfasta_file, exon_file):
    
    cmd = f"trainGlimmerHMM {mfasta_file} {exon_file}"
    # run command in shell
    run_command(cmd)

def get_bad_genes(etraining_err_file, bad_gene_file):
    
    fi  = open(etraining_err_file, "r")
    fo  = open(bad_gene_file, "w")
    for line in fi:
        pattern = re.compile("n sequence (\S+):.*")
        record = re.search(pattern, line).group(1)
        fo.write(record + "\n")
    fi.close()
    fo.close()

def count_gene_model(gb_file):
    with open(gb_file, "r") as fh:
        count = 0
        for line in fh:
            if "LOCUS" in line:
                count += 1
    return count

def prepare_augustus_training_data(gff3_file, genome_fasta, species_name, 
                                   gene_number = 1000,
                                   intergenic_size=1000, config_path="augustus/config"):
    
    default_config_path = os.environ.get("AUGUSTUS_CONFIG_PATH")
    if default_config_path is None:
        # throw error
        sys.exit("AUGUSTUS_CONFIG_PATH is not set")
    # create the config directory    
    os.makedirs(config_path, exist_ok=True)
    
    # copy the augustus config files to the new directory
    directories_to_copy = ["cgp", "extrinsic", "model", "profile", "parameters"]
    for directory in directories_to_copy:
        src_dir_path = os.path.join(default_config_path, directory)
        dest_dir_path = os.path.join(config_path, directory)
        if os.path.exists(src_dir_path) and os.path.isdir(src_dir_path):
            # Copy the directory
            shutil.copytree(src_dir_path, dest_dir_path, dirs_exist_ok=True)  # dirs_exist_ok=True allows overwriting
        else:
            print(f"Source directory does not exist or is not a directory: {src_dir_path}")
    
    # create species directory for generic files
    species_dir = os.path.join(config_path, "species")
    os.makedirs(species_dir, exist_ok=True)
    src_species_dir = os.path.join(default_config_path, "species", "generic")
    dest_species_dir = os.path.join(species_dir, "generic")
    if os.path.exists(src_species_dir) and os.path.isdir(src_species_dir):
        shutil.copytree(src_species_dir, dest_species_dir, dirs_exist_ok=True)
    else:
        print(f"Source directory does not exist or is not a directory: {src_species_dir}")

    # create species directory for new species
    cmd = f"new_species.pl --species={species_name} --AUGUSTUS_CONFIG_PATH={config_path}"
    run_command(cmd)
    
    # create the training data
    genbank_file = "training.gb"
    cmd = f"gff2gbSmallDNA.pl {gff3_file} {genome_fasta} {intergenic_size} {genbank_file}"
    run_command(cmd)
    
    # run etraining
    cmd = f"AUGUSTUS_CONFIG_PATH={config_path} ;  etraining --species={species_name} {genbank_file} "
    run_command(cmd, stdout_file="etraining.log", stderr_file="etraining.stderr")

    # filter out bad genes
    bad_gene_file = "bad_gene.lst"
    get_bad_genes("etraining.stderr", bad_gene_file)
    filter_gene_genbank_file = "training.f.gb"
    cmd = f"filterGenes.pl {bad_gene_file} {genbank_file}"
    run_command(cmd, stdout_file=filter_gene_genbank_file)

    # random split the gene model
    cmd = f"randomSplit.pl {filter_gene_genbank_file} {gene_number}"
    run_command(cmd)
    
    return filter_gene_genbank_file + ".test"

def train_augustus_model(species_name, train_gb, config_path, intergenic_size=1000, threads=8):
    config_path = os.path.realpath(config_path)
    cmd = f"export AUGUSTUS_CONFIG_PATH={config_path};\
        autoAugTrain.pl -v -v -v --trainingset={train_gb} --species={species_name} --flanking_DNA={intergenic_size} --useexisting --optrounds=1 --cpus={threads}"
    run_command(cmd, stdout_file="autoAugTrain.log", stderr_file="autoAugTrain.stderr")


if __name__ == '__main__':
    """
    mamba install -c bioconda snap glimmerhmm augustus
    """
    import argparse
    parser = argparse.ArgumentParser(description="Train gene prediction model")
    parser.add_argument("model", type=str, choices=["snap", "glimmer", "augustus"], help="gene prediction model")
    parser.add_argument("gff3", type=str, help="gff3 file")
    parser.add_argument("genome", type=str, help="genome fasta file")
    parser.add_argument("--prefix", type=str, help="output model prefix", default="model")
    parser.add_argument("--gene_number", type=int, help="number of genes to use for training", default=1000)
    parser.add_argument("--intergenic_size", type=int, help="intergenic size", default=1000)

    args = parser.parse_args()

    if args.model == "snap":
        # check snap is excutable
        train_snap_model(args.gff3, args.genome, args.prefix, args.intergenic_size)
    elif args.model == "glimmer":
        prepare_glimmer_training_data(args.gff3, args.genome, args.gene_number, args.prefix)
        train_glimmer_model(f"{args.prefix}.mfasta", f"{args.prefix}.exon_coords")
    elif args.model == "augustus":
        test_gb = prepare_augustus_training_data(args.gff3, args.genome, args.prefix, args.gene_number, args.intergenic_size)
        train_augustus_model(args.prefix, test_gb, "augustus/config", args.intergenic_size)
    else:
        print("Unknown model")
        sys.exit(1)