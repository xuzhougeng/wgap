# assembly the RNA-seq data
import os

# prepare the input file 
def get_ngs_files_from_sampleTable(wildcards):
    sample = wildcards.sample
    #print(sampleTable)
    samples = sampleTable.loc[ sample, ['fq1', 'fq2']  ].to_list()
    #print(samples)
    return samples

def get_sample_of_ngs(wildcards):
    tissue = wildcards.tissue
    samples = sampleTable.loc[(sampleTable['tissue'] == tissue) & (sampleTable['technology'] == "ngs" ) ].index.to_list()
    return [ "rna-seq/ngs/03-gtf-assembly/{sample}_ngs.gtf".format(sample = sample ) for sample in samples]

def get_sample_of_ngs_ss_rf(wildcards):
    tissue = wildcards.tissue
    samples = sampleTable.loc[(sampleTable['tissue'] == tissue) & (sampleTable['technology'] == "ngs-rf" )].index.to_list()
    #print([ "rna-seq/ngs/03-gtf-assembly/{sample}_ngs-rf.gtf".format(sample = sample ) for sample in samples])
    return [ "rna-seq/ngs/03-gtf-assembly/{sample}_ngs-rf.gtf".format(sample = sample ) for sample in samples]


def get_sample_of_ngs_ss_fr(wildcards):
    tissue = wildcards.tissue
    samples = sampleTable.loc[(sampleTable['tissue'] == tissue) & (sampleTable['technology'] == "ngs-fr" )].index.to_list()
    return [ "rna-seq/ngs/03-gtf-assembly/{sample}_ngs-fr.gtf".format(sample = sample ) for sample in samples]


if SKIP_QC:
    rule pre_processing:
       input: get_ngs_files_from_sampleTable
       output:
           r1 = temp("rna-seq/ngs/01-clean-data/{sample}_R1.fq.gz"),
           r2 = temp("rna-seq/ngs/01-clean-data/{sample}_R2.fq.gz")
       run:
           import os
           os.symlink(input[0], output.r1)
           os.symlink(input[1], output.r2)
else:
    rule pre_processing:
        input: get_ngs_files_from_sampleTable
        output:
            r1 = temp("rna-seq/ngs/01-clean-data/{sample}_R1.fq.gz"),
            r2 = temp("rna-seq/ngs/01-clean-data/{sample}_R2.fq.gz")
        params:
            json = lambda wildcards : os.path.join('qc', wildcards.sample + '.json'),
            html = lambda wildcards : os.path.join('qc', wildcards.sample + '.html')
        threads: 8
        shell:
            """fastp -w {threads} -i {input[0]} -I {input[1]} \
            -o {output.r1} -O {output.r2} -j {params.json} -h {params.html}"""

# build STAR index
rule build_STAR_index:
    input: reference
    output: directory('rna-seq/ngs/star_index')
    params:
        opts = config.get("star_idx_opts", "")
    threads: config.get('star_idx_threads', 10)
    shell:"""
    STAR \
    --runThreadN {threads} \
    --runMode genomeGenerate \
    --genomeDir {output} \
    {params.opts} \
    --genomeFastaFiles {input}
    """

# align the fastq with STAR aligner
rule rna_seq_align:
    input:
        r1 = "rna-seq/ngs/01-clean-data/{sample}_R1.fq.gz",
        r2 = "rna-seq/ngs/01-clean-data/{sample}_R2.fq.gz",
        index = 'rna-seq/ngs/star_index'
    #output: temp("rna-seq/ngs/02-read-align/{sample}_Aligned.sortedByCoord.out.bam")
    output: "rna-seq/ngs/02-read-align/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        prefix = lambda wildcards: os.path.join("rna-seq/ngs/02-read-align", wildcards.sample + "_")
    threads: 20
    shell:"""
    STAR \
        --genomeDir {input.index} \
        --runThreadN {threads} \
        --readFilesIn {input.r1} {input.r2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN 10 \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical
    """

# assembly each tissue with stringtie
# normal sequencing
rule ngs_stringtie_assembly:
    input: "rna-seq/ngs/02-read-align/{sample}_Aligned.sortedByCoord.out.bam"
    output:"rna-seq/ngs/03-gtf-assembly/{sample}_ngs.gtf"
    params:
        opts = config.get("ngs_stringtie_assembly_opts", "")
    threads: config.get('stringtie_assembly_opts', 10)
    shell: "stringtie {params.opts} -p {threads} {input} -o {output}"

rule ngs_stringtie_merge:
    input: get_sample_of_ngs
    output: "rna-seq/ngs/04-final/{tissue}_ngs.gtf"
    params:
        label = lambda wildcards: "{tissue}".format(tissue=wildcards.tissue)
    shell: "stringtie -l {params.label} --merge -o {output} {input}"

# strand-specific seqeuncing, assume stranded library rf-firststrand
rule ngs_ss_rf_stringtie_assembly:
    input: "rna-seq/ngs/02-read-align/{sample}_Aligned.sortedByCoord.out.bam"
    output:"rna-seq/ngs/03-gtf-assembly/{sample}_ngs-rf.gtf"
    params:
        opts = config.get("ngs_fr_stringtie_assembly_opts", "")
    threads: config.get('stringtie_thread', 10)
    shell: "stringtie {params.opts} --rf -p {threads} {input} -o {output}"

rule ngs_ss_rf_stringtie_merge:
    input: get_sample_of_ngs_ss_rf
    output: "rna-seq/ngs/04-final/{tissue}_ngs-rf.gtf"
    params:
        label = lambda wildcards: "{tissue}_rf".format(tissue=wildcards.tissue)
    shell: "stringtie -l {params.label} --merge -o {output} {input}"

# strand-specific seqeuncing, assume stranded library fr-secondstrand
rule ngs_ss_fr_stringtie_assembly:
    input: "rna-seq/ngs/02-read-align/{sample}_Aligned.sortedByCoord.out.bam"
    output:"rna-seq/ngs/03-gtf-assembly/{sample}_ngs-fr.gtf"
    params:
        opts = config.get("ngs_rf_stringtie_assembly_opts", "")
    threads: config.get('stringtie_thread', 10)
    shell: "stringtie {params.opts} --fr -p {threads} {input} -o {output}"

rule ngs_ss_fr_stringtie_merge:
    input: get_sample_of_ngs_ss_fr
    output: "rna-seq/ngs/04-final/{tissue}_ngs-fr.gtf"
    params:
        label = lambda wildcards: "{tissue}_fr".format(tissue=wildcards.tissue)
    shell: "stringtie -l {params.label} --merge -o {output} {input}"

# generate the output
rule prepare_ngs_maker_gff_input:
    input: "rna-seq/ngs/04-final/{tissue}_{tech}.gtf"
    output: "rna-seq/ngs/04-final/{tissue}_{tech}.gff3"
    shell: """
    gffread -E {input} -o- | sed "s#transcript#match#g" |sed "s#exon#match_part#g" > {output}
    """


