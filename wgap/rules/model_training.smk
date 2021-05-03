
# model training module
include: "train_augustus_model.smk"
include: "train_snap_model.smk"


# snap_rounds, augustus_rounds
def get_gene_models(wildcards):
    round = wildcards.get('round', '1')
    gene_models = []
    if round in snap_rounds:
        gene_models.append( "status/snap_train_round{round}.done".format(round=round) )
    if round in augustus_rounds:
        gene_models.append( "status/augustus_train_round{round}.done".format(round=round) ) 
    return gene_models

rule model_training_base:
    input: get_gene_models
    output: touch("status/model_training_base.done")    

rule model_training_loop:
    input: get_gene_models
    output: touch("status/model_training_round{round}.done")

