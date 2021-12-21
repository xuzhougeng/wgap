

def get_augustus_config_path(wildcards):
    import os
    import shutil
    import sys
    augustus_path = shutil.which('augustus')
    if not os.access(augustus_path, os.X_OK):
       print("{} is not exectuable".format(augustus_path), file=sys.stderr)
       sys.exit(1)
    
    default_config_path = os.path.join( os.path.dirname(os.path.dirname(augustus_path)), "config")

    augustus_config_path = os.environ.get('AUGUSTUS_CONFIG_PATH', default_config_path)
    augustus_species_path = os.path.join(augustus_config_path, 'species')
    
    if not os.access(augustus_species_path, os.W_OK):
       print("{} is not exectuable".format(augustus_species_path), file=sys.stderr)
       sys.exit(1)

    return os.path.abspath(augustus_config_path)

def get_training_config_path(wildcards):
    import os
    rel_path = "gene_model/augustus/round{}/config".format(wildcards.round)
    return os.path.abspath(rel_path)

# compute teh flank region size of maker gff
def compute_flank_region_size(wildcards):
    import re
    import os 
    import math

    round = int(wildcards.round)
    gff = "gene_model/augustus/round{round}/maker.gff".format(round=round)
    if not os.path.exists(gff):
        return 

    genes = dict()
        
    for line in open(gff, "r"):
        if line.startswith("#"):
            continue
        line = line.strip()    
        gtf_line = line.split('\t')
        if len(gtf_line) != 9:
            continue
        if gtf_line[2] == 'CDS' :
            pattern = re.compile(';Parent=([^;]+)')
            gene = re.search(pattern, gtf_line[8]).group(1)
            if gene not in genes:
                genes[gene] = {}
            
            min_pos = min( int(gtf_line[3]), int(gtf_line[4]) )
            max_pos = max( int(gtf_line[3]), int(gtf_line[4]) )

            if 'start' not in genes[gene]:
                genes[gene]['start'] = min_pos
            elif genes[gene]['start'] >  min_pos:
                genes[gene]['start'] = min_pos
            
            if 'end' not in genes[gene]:
                genes[gene]['end'] = max_pos
            elif genes[gene]['end'] < max_pos:
                genes[gene]['end'] = max_pos

    # compute the average length of gene
    nGenes = 0
    totalLen = 0
    avgLen = 0

    for key in genes.keys():
        nGenes += 1
        totalLen += (genes[key]['end'] - genes[key]['start'] + 1)

    avgLen = totalLen / nGenes
    flank_size = min( math.floor(avgLen / 2 ), 10000 )
    if flank_size < 0:
        print("""
        #*********
        # WARNING: flanking_DNA has the value {} , which is smaller than 0. 
        Something must have gone wrong, there. 
        Replacing by value 10000.
        #*********""".format(flank_size), file = sys.stderr)
        flank_size = 10000

    return flank_size

# get input gff 
def get_augustus_train_input(wildcards):
    round = int(wildcards.round)

    if round == 1:
        return "maker/maker_base.gff"
    elif round > 1:
        return "maker/maker_round{round}.gff".format(round=round)
    else:
        raise ValueError("loop numbers must be 1 or greater: received %s" % wildcards.round)

def get_augustus_train_input_protein(wildcards):
    round = int(wildcards.round)

    if round == 1:
        return "maker/maker_base.all.maker.proteins.fasta"
    elif round > 1:
        return "maker/maker_round{round}.all.maker.proteins.fasta".format(round=round)
    else:
        raise ValueError("loop numbers must be 1 or greater: received %s" % wildcards.round)

# create species_dir
# output: genemodel/augustus/round{round}/species/{specie_name}
rule create_species_dir:
    params:
        augustus_config_path = get_augustus_config_path,
        training_config_path = get_training_config_path,
        specie=specie_name
    output: touch("status/augustus_round{round}_specie_dir.done")
    shell:"""
    mkdir -p gene_model/augustus/round{wildcards.round}/config &&
    cp -r {params.augustus_config_path}/{{cgp,extrinsic,model,profile}} gene_model/augustus/round{wildcards.round}/config && 
    mkdir -p gene_model/augustus/round{wildcards.round}/config/species &&
    cp -r {params.augustus_config_path}/species/generic gene_model/augustus/round{wildcards.round}/config/species &&
    new_species.pl --species={params.specie} --AUGUSTUS_CONFIG_PATH={params.training_config_path}
    """

# get high quality gene model
rule get_high_quality_gff:
    input: get_augustus_train_input
    output: "gene_model/augustus/round{round}/maker.gff"
    params:
        script_dir = script_dir
    shell:"""
    python {params.script_dir}/maker_filter.py -e 1 -d 0 {input} > {output}
    """

rule init_training_set:
    input: 
        gff = "gene_model/augustus/round{round}/maker.gff",
        genome = reference
    params:
        flank_size = compute_flank_region_size
    output: temp("gene_model/augustus/round{round}/training.gb")
    shell:"""
    gff2gbSmallDNA.pl {input.gff} {input.genome} {params.flank_size} {output}
    """

# # first tranning
rule first_etraining:
    input: 
        "gene_model/augustus/round{round}/training.gb",
        "status/augustus_round{round}_specie_dir.done"
    params:
        training_config_path = lambda wildcards : "gene_model/augustus/round{}/config".format(wildcards.round),
        specie=specie_name
    output: temp("gene_model/augustus/round{round}/etraining.stderr")
    log:
        "log/augustus_round{round}_etraining.log"
    shell:"""
    export AUGUSTUS_CONFIG_PATH={params.training_config_path} && 
    etraining --species={params.specie} {input[0]} 1> {log} 2> {output} 
    """

rule get_bad_gene_list:
    input: "gene_model/augustus/round{round}/etraining.stderr"
    output: temp("gene_model/augustus/round{round}/etraining.bad.lst")
    run:
        import re
        fi  = open(input[0], "r")
        fo  = open(output[0], "w")
        for line in fi:
            pattern = re.compile("n sequence (\S+):.*")
            record = re.search(pattern, line).group(1)
            fo.write(record + "\n")
        fi.close()
        fo.close()


rule filter_bad_gene:
    input: 
        "gene_model/augustus/round{round}/training.gb",
        "gene_model/augustus/round{round}/etraining.bad.lst"
    output: temp("gene_model/augustus/round{round}/training.f.gb")
    shell:"""
    filterGenes.pl {input[1]} {input[0]} 1> {output}
    """

rule create_aa_db:
    input: 
        genbank = "gene_model/augustus/round{round}/training.f.gb",
        protein = get_augustus_train_input_protein     
    output:
        temp("tmp/protein.fasta"),
        temp("tmp/protein.fasta.phr"),
        temp("tmp/protein.fasta.pin"),
        temp("tmp/protein.fasta.psq")
    shell:"""
    perl -n -e '$_ =~/\/gene=\"(\S+)\"/ ;print "$1\n"' {input.genbank} |
        seqkit grep -f - {input.protein}   > {output[0]} &&  
        makeblastdb -in {output[0]} -dbtype prot 
    """

rule self_vs_selfblastp:
    input:
        "tmp/protein.fasta"
    output:
        temp("tmp/blastp.txt")
    threads: 20
    shell:"""
    blastp -outfmt 6 -query {input} -db {input} -num_threads {threads} > {output}
    """


rule get_redudant_gene:
    input: "tmp/blastp.txt"
    output: temp("gene_model/augustus/round{round}/redudant.gene.lst")
    params:
        maxid = 0.7
    run:
        import pandas as pd
        import numpy as np
        df = pd.read_csv(input[0],sep="\t", header=None)
        df = df.loc[df[0] != df[1], :]
        df = df.loc[df[2] > params['maxid'],: ]
        gene_list = np.unique(df[0])
        with open(ouput[0], "w") as f:
            for gene in gene_list:
                f.writelines(gene+"\n")

rule remove_redudant_gene:
    input:
        "gene_model/augustus/round{round}/training.f.gb",
        "gene_model/augustus/round{round}/redudant.gene.lst"
    output: temp("gene_model/augustus/round{round}/training.ff.gb")
    shell:"""
    filterGenes.pl {input[1]} {input[0]} > {output} 
    """


rule downsample_gene:
    input: "gene_model/augustus/round{round}/training.ff.gb"
    output: temp("gene_model/augustus/round{round}/training.fff.gb")
    params:
        threshold = config.get("training_augustus_gene_num", 1000)
    shell:"""
    if [[ {params.threshold} -eq 0 ]]; then
        cp {input} {output}
    else
        randomSplit.pl {input} {params.threshold} &&
        mv {input}.test {output} && rm -f {input}.train
    fi
    """


rule get_final_train_geneset:
    input: "gene_model/augustus/round{round}/training.fff.gb"
    output: "gene_model/augustus/round{round}/train.gb"
    shell:"""
    cp {input} {output}
    """

rule auto_training:
    input: "gene_model/augustus/round{round}/train.gb"
    output: "gene_model/augustus/round{round}/autoAugTrain.log"
    params:
        training_config_path = get_training_config_path,
        specie=specie_name,
        optround=config.get("training_augustus_opt_round",1) ,
        flank_size=compute_flank_region_size
    shell:"""
    export AUGUSTUS_CONFIG_PATH={params.training_config_path} && 
    autoAugTrain.pl -v -v -v --trainingset={input} --species={params.specie} \
      --flanking_DNA={params.flank_size} \
      --useexisting \
      --optrounds={params.optround} 1> {output}
    """

rule augustus_model_train_status:
    input: "gene_model/augustus/round{round}/autoAugTrain.log"
    params:
        augustus_config_path = get_augustus_config_path,
        training_config_path = get_training_config_path,
        specie=specie_name
    output: touch("status/augustus_train_round{round}.done")
    run:
        import shutil
        import os
        src_dir = os.path.join(params.training_config_path, 'species', str(params.specie)) 
        des_dir = os.path.join(params.augustus_config_path, 'species', str(params.specie))
        if os.path.exists( des_dir ):
            shutil.rmtree( des_dir )
        shutil.copytree( src_dir, des_dir )
