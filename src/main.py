import os
import sys
import logging
import click
from pathlib import Path

from wgap import __version__
from wgap import tools

# set the log basic config 
logging.basicConfig(
    level = logging.INFO,
    datefmt="%Y-%m-d %H:%M",
    format="[%(asctime)s %(levelname)s] %(message)s"
)

# define exception message function
def log_exception(msg):
    logging.critical(msg)
    #logging.info("Documentation is availabe at ")
    logging.info("Issue can be raised at: https://github.com/xuzhuogeng/wgap/issues")
    sys.exit(1)

def log_info(msg):
    logging.info(msg)

## define command line
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """WGAP - a whole genome annotation pipeline
    """

#####################################################
########### model training command ##################
#####################################################    
from wgap.model_prepare import prepare_training_data
from wgap.model_train import prepare_snap_training_data, train_snap_model
from wgap.model_train import prepare_augustus_training_data, train_augustus_model
from wgap.model_train import prepare_glimmer_training_data, train_glimmer_model

@cli.group()
def train():
    """Train gene prediction models"""


# Define the snap training subcommand
@train.command('snap')
@click.argument('gff3', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.option('--prefix', default='model', help='Output model prefix')
@click.option('--intergenic_size', default=1000, type=int, help='Intergenic size')
def train_snap(gff3, genome, prefix, intergenic_size):
    """Train SNAP model"""
    zff_file = prepare_snap_training_data(genome, gff3, prefix)
    train_snap_model(zff_file, genome, prefix, intergenic_size)
    
# Define the glimmer training subcommand
@train.command('glimmer')
@click.argument('gff3', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.option('--prefix', default='model', help='Output model prefix')
@click.option('--gene_number', default=1000, type=int, help='Number of genes to use for training')
def train_glimmer(gff3, genome, prefix, gene_number):
    """Train Glimmer model"""
    prepare_glimmer_training_data(gff3, genome, gene_number, prefix)
    train_glimmer_model(f"{prefix}.mfasta", f"{prefix}.exon_coords", prefix)

# Define the augustus training subcommand
@train.command('augustus')
@click.argument('gff3', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.option('--prefix', default='model', help='Output model prefix')
@click.option('--gene_number', default=1000, type=int, help='Number of genes to use for training')
@click.option('--intergenic_size', default=1000, type=int, help='Intergenic size')
def train_augustus(gff3, genome, prefix, gene_number, intergenic_size):
    """Train AUGUSTUS model"""
    # Your AUGUSTUS training logic here
    test_gb = prepare_augustus_training_data(gff3, genome, prefix, gene_number, intergenic_size)
    train_augustus_model(prefix, test_gb, "augustus/config", intergenic_size)


#####################################################
########### model prediction command ################
#####################################################  
@cli.group()
def predict():
    """Subcommand for prediction tools"""

from wgap.predict import run_all_augustus, run_all_snap, run_all_glimmer
from wgap.predict import merge_augustus_files, merge_snap_files, merge_glimmer_files

@predict.command('augustus')
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False))
@click.argument('output_dir', type=click.Path())
@click.option('--AUGUSTUS_CONFIG_PATH', required=True, help="AUGUSTUS config path")
@click.option('--species', required=True, help="Species name for AUGUSTUS")
def augustus(input_dir, output_dir, AUGUSTUS_CONFIG_PATH, species):
    """Run AUGUSTUS on fasta files"""
    """Run AUGUSTUS on fasta files"""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    fasta_files = list(input_dir.glob('*.fa'))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    run_all_augustus(fasta_files, AUGUSTUS_CONFIG_PATH, species, output_dir)
    merge_augustus_files(output_dir, "augustus.gff3")

@predict.command('snap')
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False))
@click.argument('output_dir', type=click.Path())
@click.option('--snap_hmm_file', required=True, help="HMM file for SNAP")
def snap(input_dir, output_dir, snap_hmm_file):
    """Run SNAP on fasta files"""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    fasta_files = list(input_dir.glob('*.fa'))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    run_all_snap(fasta_files, snap_hmm_file, output_dir)
    merge_snap_files(output_dir, "snap.gff3")

@predict.command('glimmer')
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False))
@click.argument('output_dir', type=click.Path())
@click.option('--glimmerhmm_model_dir', required=True, help="Model directory for GLIMMER")
def glimmer(input_dir, output_dir, glimmerhmm_model_dir):
    """Run GlimmerHMM on fasta files"""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    fasta_files = list(input_dir.glob('*.fa'))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    run_all_glimmer(fasta_files, glimmerhmm_model_dir, output_dir)
    merge_glimmer_files(output_dir, "glimmer.gff3")


# add utr with new click
from wgap.utr import add_utr_to_gene_model
@cli.command('add_utr')
@click.argument('gff3', type=click.Path(exists=True))
#gtf_files: a list of gtf files
@click.argument('gtf_files', nargs=-1, type=click.Path(exists=True))
@click.option('--output', default='output.gff3', help='Output file name')
def add_utr(gff3, gtf_files, output):
    """Add UTR to GFF3 file"""
    # Your add UTR logic here
    update_genes = add_utr_to_gene_model(gff3, gtf_files)
    with open(output, 'w') as f:
        for gene in update_genes:
            f.write(gene.to_gff3())




#####################################################
########### other useful tool command ###############
#####################################################

@cli.group()
def tool():
    """Subcommand for prediction tools"""

# download command for homology evidence
@tool.command(
    'download',
    context_settings = {"ignore_unknown_options" : True},
    short_help = "dowanlod the protein sequence from Swiss/UniPort"
)
@click.argument(
    "fasta",
    default="protein.fa"
)
@click.option(
    "-s", "--specie",
    default="plants",
    type = click.Choice(["archaea", "bacteria", "fungi", "human", "invertebrates",
    "mammals","plants","rodents", "vertebrates", "viruses"]),
    help = "choose your own specie"
)
@click.option(
    "-d", "--dataset",
    default="sprot",
    type=click.Choice(["sprot", "trembl"]),
    help="sprot is recommended"
)
def download_protein(fasta, specie, dataset):

    logging.info("Downloading: %s" % dataset + "_"+ specie )
    dat_file = ""
    dat_file = tools.download_uniprot(specie, dataset)
    if not dat_file:
        sys.exit(1)
    tools.convet_dat_to_fasta(dat_file, fasta)
    logging.info("Finished: %s" % fasta )


if __name__ == '__main__':
    cli()
