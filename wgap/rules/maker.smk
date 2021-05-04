
# rule to generate maker.gff
rule maker_annotation_output:
    input: "status/maker_annotation.done"
    params:
         base_name = "maker",
         dsindex = "{ds}.maker.output/{ds}_master_datastore_index.log".format(ds=specie_name)
    output: "maker.gff"
    shell:"""
    fasta_merge -o {params.base_name} -d {params.dsindex} && 
    gff3_merge -o {output} -d {params.dsindex}
    """

# CTL file rules
def add_param_to_setting(settings, params):
    if len(params) == 0:
        return settings
    for param in params:
        key,value = param.split("=")
        settings[key] = value
    return settings

def get_maker_config(wildcards):

    maker_ctl  = wildcards.ctl
    maker_type = wildcards.type
    json_file = ".maker_{}_{}_parameter.json".format(maker_ctl, maker_type)
    settings = dict()

    # maker round
    if maker_type.find("round") >= 0:
        round = int(maker_type[5:])
    elif maker_type == "base":
        round = 1
    else:
        round = -1 # prediction
    # get the params
    if round > 0:
        params = config['round'+str(round)].get(maker_ctl,"")
    else:
        params = config['prediciton'].get(maker_ctl,"")
        
    if maker_ctl == "bopts" or maker_ctl == "exe": # bopts and exe
        settings =  add_param_to_setting(settings, params)
    else: # opts setting
        # default settings
        settings['genome'] = reference
        if transcript_assemble:
            settings['est_gff'] = ','.join(all_ngs_gff3)
        else:
            settings['est_gff'] = gtf
        settings['protein'] = protein
        if maker_type == "annotation":
            # training model: setting snap, hmm 
            if training_model:
                if snap_rounds != '0':
                    settings['snaphmm'] = "gene_model/snap_round{}.hmm".format(snap_rounds[-1])
                if augustus_rounds != '0':
                    settings["augustus_species"] = specie_name
            else: 
                settings['snaphmm'] = snaphmm
                settings['augustus_species'] = augustus_species


        if maker_type == "base": # train round = 1
            if snaphmm or augustus_species:
                settings['snaphmm'] = snaphmm
                settings['augustus_species'] = augustus_species
            else:
                settings['est2genome'] = '1'
                settings['protein2genome'] = '1'
        
        if maker_type.find("round") >= 0: #  train round > 1
            if str(round-1) in snap_rounds :
                settings['snaphmm'] = "gene_model/snap_round{}.hmm".format(round-1)
            if str(round-1) in augustus_rounds :
                    settings["augustus_species"] = specie_name

        settings =  add_param_to_setting(settings, params)
    # dump the settings to json
    with open(json_file, "w") as f:
        json.dump(settings, f)
    return os.path.abspath(json_file)


# return the input base on the training_model and transcript_assemble
def get_maker_ctl_input(wildcards):

    maker_ctl  = wildcards.ctl
    maker_type = wildcards.type
    
    input_dict = {
        "ctl" : "maker_{ctl}.ctl".format(ctl=maker_ctl)
    }

    if maker_ctl == "bopts" or maker_ctl == "exe":
        return input_dict

    if maker_type == "annotation":
        if training_model:
            input_dict['train'] = "status/model_training_round{round}.done".format(round=max_round)
            #print(input_dict)
            return input_dict
        elif transcript_assemble:
            input_dict['gtf'] = "status/transcript_assemble.done"
            #print(input_dict)
            return input_dict
        elif check_gff(gtf):
            return input_dict
        else:
            log_exception("gtf not exists")

    if maker_type == "base":
        #print(input_dict)
        if transcript_assemble:
            input_dict['gtf'] = "status/transcript_assemble.done"
            return input_dict
        elif check_gff(gtf):
            return input_dict
        else:
            log_exception("gtf not exists")

    if maker_type.find("round") >= 0:
        round = int(maker_type[5:])
        input_dict['train'] = "status/model_training_round{round}.done".format(round=round-1)
        #print(input_dict)
        return input_dict

rule init_maker_ctl:
    output: "maker_opts.ctl", "maker_bopts.ctl", "maker_exe.ctl"
    shell:"maker -CTL"

rule update_maker_ctl:
    input: unpack(get_maker_ctl_input)
    params:
        script_dir = script_dir,
        opts_settings = get_maker_config
    output: "ctl/maker_{ctl}_{type}.ctl"
    shell:"""
    python {params.script_dir}/maker_ctl_updator.py {input.ctl} {output} {params.opts_settings}
    """

# annotation step
if MPI:
    rule run_maker_annotation:
        input: 
            opts = "ctl/maker_opts_annotation.ctl",
            bopts = "ctl/maker_bopts_annotation.ctl",
            exe = "ctl/maker_exe_annotation.ctl"
        output: touch("status/maker_annotation.done")
        params:
            dsname = specie_name
        threads: 40
        log: "log/maker_annotation.log"
        shell:"""
        mpiexec -n {threads} maker -quiet -base {params.dsname} {input.opts} {input.bopts} {input.exe} &> {log} 
        """
else: 
    rule run_maker_annotation:
        input: 
            opts = "ctl/maker_opts_annotation.ctl",
            bopts = "ctl/maker_bopts_annotation.ctl",
            exe = "ctl/maker_exe_annotation.ctl"
        output: touch("status/maker_annotation.done")
        params:
            dsname = specie_name,
        threads: 1
        log: "log/maker_annotation.log"
        shell:"""
        maker -quiet -base {params.dsname} {input.opts} {input.bopts} {input.exe} &> {log} 
        """

# train step: maker_loop: maker_round{2...n}.gff
rule maker_loop_output:
    input: "status/maker_round{round}.done"
    params:
        base_name =  lambda wildcards: "maker/maker_round{round}".format(round=wildcards.round),
        dsindex = "{ds}.maker.output/{ds}_master_datastore_index.log".format(ds=specie_name)
    output: 
        gff = "maker/maker_round{round}.gff",
        protein = "maker/maker_round{round}.all.maker.proteins.fasta",
        transcript = "maker/maker_round{round}.all.maker.transcripts.fasta"
    shell:"""
    fasta_merge -o {params.base_name} -d {params.dsindex} && 
    gff3_merge -o {output.gff} -d  {params.dsindex}     
    """

if MPI:
    rule run_maker_loop:
        input: 
            opts = "ctl/maker_opts_round{round}.ctl",
            bopts = "ctl/maker_bopts_round{round}.ctl",
            exe = "ctl/maker_exe_round{round}.ctl"
        output: touch("status/maker_round{round}.done")
        params:
            dsname = specie_name,
        threads: 40
        log: "log/maker_round{round}.log"
        shell:"""
        mpiexec -n {threads} maker -quiet -base {params.dsname} {input.opts} {input.bopts} {input.exe} &> {log} 
        """
else:
    rule run_maker_loop:
        input: 
            opts = "ctl/maker_opts_round{round}.ctl",
            bopts = "ctl/maker_bopts_round{round}.ctl",
            exe = "ctl/maker_exe_round{round}.ctl"
        output: touch("status/maker_round{round}.done")
        params:
            dsname = specie_name,
        threads: 1
        log: "log/maker_round{round}.log"
        shell:"""
        {threads} maker -quiet -base {params.dsname} {input.opts} {input.bopts} {input.exe} &> {log} 
        """

# maker_base.gff
if MPI:
    rule run_maker_base:
        input: 
            opts = "ctl/maker_opts_base.ctl",
            bopts = "ctl/maker_bopts_base.ctl",
            exe = "ctl/maker_exe_base.ctl"
        output: touch("status/maker_base.done")
        params:
            dsname = specie_name,
        threads: 40
        log: "log/maker_base.log"
        shell:"""
        mpiexec -n {threads} maker -quiet -base {params.dsname} {input.opts} {input.bopts} {input.exe} &> {log} 
        """
else:
    rule run_maker_base:
        input: 
            opts = "ctl/maker_opts_base.ctl",
            bopts = "ctl/maker_bopts_base.ctl",
            exe = "ctl/maker_exe_base.ctl"
        output: touch("status/maker_base.done")
        params:
            dsname = specie_name,
        threads: 1
        log: "log/maker_base.log"
        shell:"""
        maker -quiet -base {params.dsname} {input.opts} {input.bopts} {input.exe} &> {log} 
        """

rule maker_base_output:
    input: "status/maker_base.done"
    params:
        base_name = "maker/maker_base",
        dsindex = "{ds}.maker.output/{ds}_master_datastore_index.log".format(ds=specie_name)
    output: 
        gff = "maker/maker_base.gff",
        protein = "maker/maker_base.all.maker.proteins.fasta",
        transcript = "maker/maker_base.all.maker.transcripts.fasta"
    shell:"""
    fasta_merge -o {params.base_name} -d {params.dsindex} &&
    gff3_merge -o {output.gff} -d  {params.dsindex}     
    """