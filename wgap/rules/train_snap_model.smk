# snap gene model training module

def get_snap_train_input(wildcards):
    round = int(wildcards.round)

    if round == 1:
        return "maker/maker_base.gff"
    elif round > 1:
        return "maker/maker_round{round}.gff".format(round=round)
    else:
        raise ValueError("loop numbers must be 1 or greater: received %s" % wildcards.round)

rule extract_gff_from_maker:
    input: get_snap_train_input
    params:
        ss = "0.5",
        est_ovl = "0.5",
        AED = "0.5"
    output: 
        "gene_model/snap/round{round}/genome.ann", 
        "gene_model/snap/round{round}/genome.dna", 
    shell:"""
    mkdir -p  gene_model/snap/round{wildcards.round} && \
    cd gene_model/snap/round{wildcards.round} && \
    maker2zff -c {params.ss} -e {params.est_ovl} -x {params.AED} ../../../{input}
    """

rule snap_fathom_categorize:
    input:
        "gene_model/snap/round{round}/genome.ann", 
        "gene_model/snap/round{round}/genome.dna"
    output:
        "gene_model/snap/round{round}/uni.ann",
        "gene_model/snap/round{round}/uni.dna"
    shell:"""
    cd gene_model/snap/round{wildcards.round} && \
    fathom -categorize 1000 genome.ann genome.dna
    """

rule snap_fathom_export:
    input:
        "gene_model/snap/round{round}/uni.ann",
        "gene_model/snap/round{round}/uni.dna"   
    output:
        "gene_model/snap/round{round}/export.ann",
        "gene_model/snap/round{round}/export.dna"   
    shell:"""
    cd gene_model/snap/round{wildcards.round} && \
    fathom -export 1000 -plus uni.ann uni.dna 
    """

rule snap_fathom_forge:
    input:
        "gene_model/snap/round{round}/export.ann",
        "gene_model/snap/round{round}/export.dna"
    output: directory("gene_model/snap/round{round}/params")   
    shell:"""
    cd gene_model/snap/round{wildcards.round} && \
    mkdir -p params && cd params && \
    forge ../export.ann ../export.dna 
    """

rule snap_hmm_assemble:
    input: directory("gene_model/snap/round{round}/params")   
    params:
        species = specie_name
    output: "gene_model/snap_round{round}.hmm"
    shell:"""
    hmm-assembler.pl {params.species} {input} > {output}
    """

rule snap_model_train_status:
    input: "gene_model/snap_round{round}.hmm"
    output: touch("status/snap_train_round{round}.done")