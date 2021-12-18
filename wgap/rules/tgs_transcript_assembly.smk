# prepare the input file 
def get_tgs_files_from_sampleTable(wildcards):
    sample = wildcards.sample
    #print(sampleTable)
    sample = sampleTable.loc[ sample, 'fq1' ]
    #print(samples)
    return sample

def get_sample_of_pb(wildcards):
    tissue = wildcards.tissue
    samples = sampleTable.loc[(sampleTable['tissue'] == tissue) & (sampleTable['technology'] == "pb" ) ].index.to_list()
    return [ "rna-seq/tgs/03-gtf-assembly/{sample}_pb.gtf".format(sample = sample ) for sample in samples]

def get_sample_of_ont(wildcards):
    tissue = wildcards.tissue
    samples = sampleTable.loc[(sampleTable['tissue'] == tissue) & (sampleTable['technology'] == "ont" ) ].index.to_list()
    return [ "rna-seq/tgs/03-gtf-assembly/{sample}_ont.gtf".format(sample = sample ) for sample in samples]


rule build_mm2_index:
    input: reference
    output: 'rna-seq/tgs/genome_index.mmi'
    params:
        opts = config.get("mm2_idx_opts", "")
    threads: config.get('mm2_idx_threads', 10)
    shell:"""
        minimap2 -t {threads} {params.opts} -I 100G -d {output} {input}
    """

rule tgs_pre_processing:
    input: get_tgs_files_from_sampleTable
    output: temp("rna-seq/tgs/01-clean-data/{sample}.fq.gz"),
    run:
        import os
        os.symlink(input, output)

rule pb_map:
    input:
        index = rules.build_mm2_index.output,
        fastq = rules.tgs_pre_processing.output
    output:
        bam = "rna-seq/tgs/02-read-align/{sample}_pb.bam"
    params: 
        opts = config.get("mm2_pb_map_opts", ""),
        min_mq = config.get("mm2_pb_min_mq", 40),
    threads: config.get("mm2_map_threads", 20)
    shell:"""
    minimap2 -t {threads} {params.opts} -ax splice:hq -uf {input.index} {input.fastq} \ |
      samtools view -q {x.min_mq} -F 2304 -Sb - \
      samtools sort -@ {threads} -o {output.bam} - && \
      samtools index {output.bam}  
    """


rule pb_ont:
    input:
        index = rules.build_mm2_index.output,
        fastq = rules.tgs_pre_processing.output
    output: 
        bam = "rna-seq/tgs/02-read-align/{sample}_ont.bam"
    params: 
        opts = config.get("mm2_ont_map_opts", ""),
        min_mq = config.get("mm2_ont_min_mq", 40),
    threads: config.get("mm2_map_threads", 20)
    shell:"""
    minimap2 -t {threads} {params.opts} -ax splice -uf {input.index} {input.fastq} \ |
      samtools view -q {params.min_mq} -F 2304 -Sb - \
      samtools sort -@ {threads} -o {output.bam} - && \
      samtools index {output.bam}  
    """

# long read RNA sequencing
rule pb_stringtie_assembly:
    input: "rna-seq/tgs/02-read-align/{sample}_pb.bam"
    output:"rna-seq/tgs/03-gtf-assembly/{sample}_pb.gtf"
    params:
        opts = config.get("tgs_stringtie_assembly_opts", "")
    threads: config.get('stringtie_assembly_opts', 10)
    shell: "stringtie {params.opts} -p {threads} {input} -o {output}"

rule pb_stringtie_merge:
    input: get_sample_of_pb
    output: "rna-seq/tgs/04-final/{tissue}_pb.gtf"
    params:
        label = lambda wildcards: "{tissue}_pb".format(tissue=wildcards.tissue)
    shell: "stringtie --merge -l {params.label} -o {output} {input}"


rule ont_stringtie_assembly:
    input: "rna-seq/tgs/02-read-align/{sample}_ont.bam"
    output:"rna-seq/tgs/03-gtf-assembly/{sample}_ont.gtf"
    params:
        opts = config.get("tgs_stringtie_assembly_opts", "")
    threads: config.get('stringtie_assembly_opts', 10)
    shell: "stringtie {params.opts} -p {threads} {input} -o {output}"

rule ont_stringtie_merge:
    input: get_sample_of_ont
    output: "rna-seq/tgs/04-final/{tissue}_ont.gtf"
    params:
        label = lambda wildcards: "{tissue}_ont".format(tissue=wildcards.tissue)
    shell: "stringtie --merge -l {params.label}  -o {output} {input}"

# generate the output
rule prepare_tgs_maker_gff_input:
    input: "rna-seq/tgs/04-final/{tissue}_{tech}.gtf"
    output: "rna-seq/tgs/04-final/{tissue}_{tech}.gff3"
    shell: """
    gffread -E {input} -o- | sed "s#transcript#match#g" |sed "s#exon#match_part#g" > {output}
    """
