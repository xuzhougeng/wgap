

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
        script_dir = script_dir,
        filter_opt = "-e 1 -d 0"
    shell:"""
    python {params.script_dir}/maker_filter.py {params.filter_opt} {input} > {output}
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

# first tranning to remove training error
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

rule self_vs_self_blastp:
    input: 
        genbank = "gene_model/augustus/round{round}/training.f.gb",
        protein = get_augustus_train_input_protein     
    output:
        orig=temp("gene_model/augustus/round{round}/training.f.aa.fasta"),
        flt=temp("gene_model/augustus/round{round}/non_reduant_gene.fasta")
    params:
        maxid = 0.8
    threads: 20
    shell:"""
    perl -n -e '$_ =~/\/gene=\"(\S+)\"/ ;print "$1\n"' {input.genbank} |
        seqkit grep -f - {input.protein} > {output.orig} &&  
        aa2nonred.pl --maxid={params.maxid} --diamond --cores={threads} {output.orig} {output.flt}
    """

rule get_nonredudant_gene:
    input: 
        aa = "gene_model/augustus/round{round}/non_reduant_gene.fasta",
        genbank = "gene_model/augustus/round{round}/training.f.gb",
    output: temp("gene_model/augustus/round{round}/training.ff.gb")      
    run:
        import re
        gene_list = set()
        with open(input['aa'], "r") as f:
            for line in f.readlines():
                if line.startswith(">"):
                    gene_list.add( line[1:] )

        buffer_line = ""
        gene_pattern = re.compile(r'\/gene=\"(\S+)\"')
        new_gb = open(output[0], "w")
        with open(input['genbank'], "r") as f:
            for line in f.readlines():
                if not line.startswith("//"):
                    buffer_line += line
                else:
                    gene_id = re.findall(gene_pattern, buffer_line)[0]
                    if gene_id not in gene_list:
                        buffer_line += "//\n"
                        new_gb.write(buffer_line)
                    buffer_line = ""
        new_gb.close()

rule downsample_gene:
    input: "gene_model/augustus/round{round}/training.ff.gb"
    output: temp("gene_model/augustus/round{round}/training.fff.gb")
    params:
        threshold = config.get("training_augustus_gene_num", 1000)
    shell:"""
    gene_num=$(grep -c "LOCUS" {input})
    if [ {params.threshold} -eq 0 ] || [ $gene_num -lt {params.threshold} ]; then
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
        flank_size=compute_flank_region_size,
        auto_aug_opts = config.get("training_augustus_opts", "")
    threads: 8
    shell:"""
    export AUGUSTUS_CONFIG_PATH={params.training_config_path} && 
    autoAugTrain.pl -v -v -v --trainingset={input} --species={params.specie} \
      --flanking_DNA={params.flank_size} \
      --useexisting \
      --optrounds={params.optround} \
      --cpus={threads} 1> {output}
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
