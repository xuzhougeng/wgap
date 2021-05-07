
from wgap.scripts.maker_update import gff_reader
from wgap.scripts.maker_update import gff_writer
from BCBio import GFF
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
#from Bio.SeqFeature import FeatureLocation


def get_model_dict(records):
    model_dict = {}
    for rec in records:
        for feature in rec.features:
            model_dict[feature.id] = feature
    return model_dict

# create gene feature
def create_gene_model(seqfeature, source):
    
    mRNA_name = seqfeature.qualifiers.get('Name')[0].replace('abinit','processed')
    gene_name = mRNA_name.replace("-mRNA-1", "")
    match = seqfeature
    match_part = match.sub_features
    # gene 
    gene_qualifiers = {
        "ID": [gene_name],
        "source": [source],
        "Name": [gene_name]
    }
    gene_feature = SeqFeature(
        match.location, type="gene", qualifiers=gene_qualifiers, id = gene_name
    )
    # mRNA
    mRNA_qualifiers = {
        "ID": [mRNA_name],
        "Parent": [gene_name],
        "Name": [mRNA_name],  
        "source": [source]
    }
    mRNA_feature = SeqFeature(
        match.location, type="mRNA", qualifiers=mRNA_qualifiers, id = mRNA_name
    )
    gene_feature.sub_features = mRNA_feature
    # CDS and exon feature
    feature_list = []
    count = 0 
    for part in match_part:
        count += 1
        exon_name = mRNA_name + ":" + str(count)
        cds_name  = mRNA_name + ":cds"
        exon_qualifiers = {
            "ID": [exon_name],
            "Parent": [mRNA_name],
            "source": [source]
        }
        cds_qualifiers = {
            "ID": [cds_name],
            "Parent": [mRNA_name],
            "source": [source],
        }
        
        cds_feature = SeqFeature(
            part.location, type="CDS", qualifiers=cds_qualifiers, id = cds_name
        )
        exon_feature = SeqFeature(
            part.location, type="exon", qualifiers=exon_qualifiers, id = exon_name
        )
        feature_list.append(cds_feature)
        feature_list.append(exon_feature)
    mRNA_feature.sub_features = feature_list

    return gene_feature



# add new gene to existed models
def add_gene_model(seqfeature, orig_models, idx):
    from bisect import bisect_left
    # format the gene model
    seqfeature = create_gene_model(seqfeature,  "maker" )
    
    # find the insert position
    position = [ f.location.nofuzzy_start for f in orig_models[idx].features]
    if len(position) == 0:
        insert_site = 0
    else:
        insert_site = bisect_left(position, seqfeature.location.nofuzzy_start)

    # logging
    #print("add gene: " + gene_name,file=stderr)
    #print(seqfeature)

    orig_models[idx].features.insert(insert_site, seqfeature)


def update_gene_model(makergff, abinitgff, geneset ):

    # get the gene and corresponding feature
    orig_models = gff_reader(makergff)
    # get reference index in SeqRecord list
    orig_index = {l.id : idx for idx,l in enumerate( orig_models )}
    max_idx = max(list(orig_index.values()))

    # read the new model from apollo exported gff
    in_handle = open(abinitgff, "r")
    for rec in GFF.parse(in_handle):
        idx = orig_index.get(rec.id,-1)
        if idx >= 0:
            for feature in rec.features:
                gene_name = feature.qualifiers.get('Name')[0]
                if gene_name not in geneset:
                    continue

                add_gene_model(feature, orig_models, idx)
        else:
            seqrecord = SeqRecord(seq=rec.seq, id=rec.id)
            orig_models.append(seqrecord) 
            max_idx += 1
            orig_index[rec.id] = max_idx
            add_gene_model(feature, orig_models, max_idx)
    
    in_handle.close()

    return orig_models


def rescue(makergff, abinitgff, genelist, outgff):

    # get the gene list
    geneset = set()

    for gene in open(genelist, "r"):
        gene = gene.strip()
        gene = gene.replace('processed', 'abinit')
        geneset.add(gene)

    out_models = update_gene_model(makergff, abinitgff, geneset)
    gff_writer(out_models, outgff)

