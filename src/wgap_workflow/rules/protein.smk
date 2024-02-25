
# exonerate --model protein2genome query.fasta target.fasta


rule miniprot_index:
    input:
        "genome.fa"
    output:
        "genome.mpi"
    shell:
        "miniprot -t 16 -d genome.mpi {input}"

rule protein2genome:
    input:
        query = "bol.faa",
        target = "genome.fa,
        index = "genome.mpi"
    output:
        "bol.pep.gff"
    shell:
        "miniprot -Iut16 --gff {input.index} {input.query} > {output}"

