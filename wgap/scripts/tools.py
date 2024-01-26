# download the protin from UniPort
import os
import sys
import gzip
import urllib.request
from wgap import root_dir
from Bio import SwissProt

def get_snakefile(file = "Snakefile"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file;  tried %s" %sf)
    return sf

def get_configfile(file = "template_config.yaml"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the config.yaml file;  tried %s" %sf)
    return sf

def get_samplefile(file = "template_sample.csv"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the sample.csv file;  tried %s" %sf)
    return sf


#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_plants.dat.gz
#zcat uniprot_sprot_plants.dat.gz |\
#    awk '{if (/^ /) {gsub(/ /, ""); print} else if (/^AC/) print ">" $2}' |\
#    sed 's/;$//'> protein.fa

def download_uniprot(specie, dataset):
    
    file_name = "uniprot_{dataset}_{specie}.dat.gz".format(dataset=dataset, specie=specie)
    file_path = file_name
    base_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/"
    url = base_url + file_name

    try:
        urllib.request.urlretrieve(url, filename=file_path)
    except:
        print("You can also download protein with `wget {}`, and then rerun wgap".format(url))
    
    return file_path


# convert SwissPort to FASTA
def convet_dat_to_fasta(dat, fasta):
    handle = gzip.open(dat, "rt")
    file_out = open(fasta, "w")
    for record in SwissProt.parse(handle):
        gene_name = record.gene_name
        organism = record.organism
        seq = record.sequence
        out_line = ">{}{}\n{}\n".format(gene_name, organism, seq)
        file_out.writelines(out_line)





