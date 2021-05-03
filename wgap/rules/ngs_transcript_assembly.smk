# assembly the RNA-seq data
import os


# prepare the input file 
def get_files_from_sampleTable(wildcards):
    sample = wildcards.sample
    #print(sampleTable)
    samples = sampleTable.loc[ sample, ['fq1', 'fq2']  ].to_list()
    #print(samples)
    return samples

def get_sample_of_same_tissue(wildcards):
    tissue = wildcards.tissue
    samples = sampleTable.loc[sampleTable['tissue'] == tissue].index.to_list()
    return [ "rna-seq/ngs/03-gtf-assembly/{sample}.gtf".format(sample = sample ) for sample in samples]


if SKIP_QC:
    rule pre_processing:
       input: get_files_from_sampleTable
       output:
           r1 = temp("rna-seq/ngs/01-clean-data/{sample}_R1.fq.gz"),
           r2 = temp("rna-seq/ngs/01-clean-data/{sample}_R2.fq.gz")
       run:
           import os
           os.symlink(input[0], output.r1)
           os.symlink(input[1], output.r2)
else:
    rule pre_processing:
        input: get_files_from_sampleTable
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
rule build_index:
    input: reference
    output: directory('rna-seq/ngs/star_index')
    params:
        SAindexNbases = "14"
    threads: 10
    shell:"""
    STAR \
    --runThreadN {threads} \
    --runMode genomeGenerate \
    --genomeDir {output} \
    --genomeSAindexNbases {params.SAindexNbases}\
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
rule stringtie_assembly:
    input: "rna-seq/ngs/02-read-align/{sample}_Aligned.sortedByCoord.out.bam"
    output:"rna-seq/ngs/03-gtf-assembly/{sample}.gtf"
    threads: 10
    shell: "stringtie -p {threads} {input} -o {output}"

rule stringtie_merge_tissue:
    input: get_sample_of_same_tissue
    output: "rna-seq/ngs/04-final/{tissue}.gtf"
    shell: "stringtie --merge -o {output} {input}"

rule prepare_maker_gff_input:
    input: "rna-seq/ngs/04-final/{tissue}.gtf"
    output: "rna-seq/ngs/04-final/{tissue}.gff3"
    shell: """
    gffread -E {input} -o- | sed "s#transcript#match#g" |sed "s#exon#match_part#g" > {output}
    """

"""
rule trinity_denovo_assembly:
    input:
    output:

rule trinity_guided_assembly:
"""


