import os
import sys
import logging
import click
from pathlib import Path
import importlib.resources as pkg_resources

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


def get_snakefile_path(file = "Snakefile"):
    """
    Snakefile, config.yaml, template_config.yaml, template_sample.csv
    """

    with pkg_resources.path('wgap_workflow', file) as snakefile_path:
        print(snakefile_path)

## define command line
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """WGAP - a whole genome annotation pipeline
    """


############################################################
###########  prepare from downstream analysis ##############
############################################################
from wgap.model_prepare import  prepare_training_data

@cli.group()
def prepare():
    """Prepare the input files for annotation"""
    

@prepare.command('transcriptome')
def transcript_assembly():
    """Prepare transcriptome assembly for gene prediction"""
    # Your transcriptome assembly logic here
    pass

@prepare.command('protein')
def prepare_protein():
    """Prepare protein sequence for gene prediction"""
    # Your protein preparation logic here
    pass

@prepare.command('gtf2train')
@click.argument('gtf', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.option('--gtf_source', default='stringtie',  help='Output gff format, default=stringtie')
@click.option('--prefix', default='wgap', help='Output model prefix')
@click.option('--min_exon_num', default=3, type=int, help='Minimum exon number, default=3')
@click.option('--max_exon_num', default=10000, type=int, help='Maximum exon number,default=10000')
@click.option('--min_orf_size', default=100, type=int, help='Minimum ORF size, default=100')
def gtf2train(gtf, genome, gtf_source, prefix, min_exon_num, max_exon_num, min_orf_size ):
    """Prepare training data for gene prediction"""
    # Your training data preparation logic here
    prepare_training_data(gtf, genome, gtf_source, prefix, min_exon_num, max_exon_num, min_orf_size )


# # workflow command
# @cli.command(
#     "run",
#     context_settings = {"ignore_unknown_options" : True},
#     short_help = "run wgap main workflow"
# )
# @click.argument(
#     "workflow",
#     default="all",
#     type=click.Choice(["all", "transcript_assembly", "gene_model_training", "test"]),
# )
# @click.option("-w",
#     "--working-dir",
#     type = click.Path(dir_okay = True, writable=True, resolve_path=True),
#     help = "project running directory",
#     default = "."
# )
# @click.option("-c",
#     "--config-file",
#     type = click.Path(dir_okay = True, writable=True, resolve_path=True),
#     help = "config file, generated with wgap init",
# )
# @click.option(
#     "-j",
#     "--jobs",
#     type=int,
#     help="use at most this many jobs in parallel (see cluster submission for mor details).",
# )
# @click.option(
#     "--profile",
#     default=None,
#     help="snakemake profile e.g. for cluster execution.",
# )
# @click.option(
#     "-n",
#     "--dryrun",
#     is_flag=True,
#     default=False,
#     show_default=True,
#     help="Test execution.",
# )
# @click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
# def run_workflow(workflow, working_dir, config_file, jobs, profile, dryrun, snakemake_args):

#     # check file availability
#     if config_file is None:
#         config_file = os.path.join(working_dir, 'config.yaml')
    
#     if not os.path.exists(config_file):
#         logging.critical(f'config-file not found: {config_file}.')
#         sys.exit(1)
    
#     conf = load_configfile(config_file)
#     sample_file = conf['sample']

#     if not os.path.exists(sample_file):
#         logging.critical('{} is not existed'.format(sample_file))
#         sys.exit(1)
    
#     # TO DO: check file avaiable
    
#     cmd = (
#         "snakemake --snakefile {snakefile} --directory {working_dir}"
#         " {jobs} --rerun-incomplete "
#         " --configfile '{config_file}' --nolock "
#         " {profile} {dryrun}"
#         " {target_rule} "
#         " {args} "
#     ).format(
#         snakefile = utils.get_snakefile(),
#         working_dir = working_dir,
#         jobs = "--jobs {}".format(jobs if jobs is not None else "4"),
#         config_file = config_file,
#         profile = "" if (profile is None) else "--profile {}".foramt(profile),
#         dryrun="--dryrun" if dryrun else "",
#         args = " ".join(snakemake_args),
#         target_rule = workflow if workflow != "None" else ""
#     )
#     logging.info("Executing: %s" % cmd)
#     try:
#         subprocess.check_call(cmd, shell=True)
#     except subprocess.CalledProcessError as e:
#         logging.critical(e)
#         exit(1)

# # initiate project
# @cli.command(
#     'init',
#     context_settings = {"ignore_unknown_options" : True},
#     short_help = "copy the config.yaml and samples.csv to working directoty"
# )
# @click.argument(
#     'workdir',
#     default="."
# )
# def init_workdir(workdir):
#     config_file = utils.get_configfile()
#     sample_file = utils.get_samplefile()

#     copyfile(config_file, os.path.join(workdir, "config.yaml"))
#     copyfile(sample_file, os.path.join(workdir, "samples.csv"))




#####################################################
########### model training command ##################
#####################################################    
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
    gff3 = os.path.abspath(gff3)
    genome = os.path.abspath(genome)
    zff_file = prepare_snap_training_data(gff3, genome, prefix)
    train_snap_model(zff_file, genome, prefix, intergenic_size)
    
# Define the glimmer training subcommand
@train.command('glimmer')
@click.argument('gff3', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.option('--prefix', default='model', help='Output model prefix')
@click.option('--gene_number', default=1000, type=int, help='Number of genes to use for training')
def train_glimmer(gff3, genome, prefix, gene_number):
    """Train Glimmer model"""
    gff3 = os.path.abspath(gff3)
    genome = os.path.abspath(genome)
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
    # get abs path
    gff3 = os.path.abspath(gff3)
    genome = os.path.abspath(genome)

    test_gb = prepare_augustus_training_data(gff3, genome, prefix, gene_number, intergenic_size)
    train_augustus_model(prefix, test_gb, "augustus/config", intergenic_size)


#####################################################
########### model prediction command ################
#####################################################  
from wgap.predict import run_all_augustus, run_all_snap, run_all_glimmer
from wgap.predict import merge_augustus_files, merge_snap_files, merge_glimmer_files
from wgap.split_genome import split_genome_by_window, split_genome_by_evidence, split_genome_by_te
from wgap.tools import add_exon_to_gene_model


@cli.group()
def predict():
    """Subcommand for prediction tools"""

@predict.command('split')
@click.argument('genome', type=click.Path(exists=True))
@click.option('--outdir', default='split_genome', help='Output directory')
@click.option('--method', default='windows', help='Split method', type=click.Choice(['windows', 'evidence', 'te']))
@click.option('--window_size', default=100000, help='Window size')
@click.option('--repeat_masker', help='Repeat masker result')
@click.option('--homology', help='Homology result in gff format')
@click.option('--gene_gff', help='Gene gff file for evidence based split')
def split_genome(genome, method, window_size, repeat_masker, gene_gff, homology, outdir):
    """Split genome into windows for parallel gene prediction"""
    # Your split genome logic here

    if method == "windows":
        split_genome_by_window(genome, window_size, outdir)
    elif method == 'evidence':
        split_genome_by_evidence(genome, gene_gff, homology, outdir)
    elif method == "te":
        split_genome_by_te(genome, repeat_masker, gene_gff, homology, 100, outdir)
    else:
        print("method not supported")
        sys.exit(1)

@predict.command('augustus')
@click.argument('input_dir', type=click.Path(exists=True, file_okay=False))
@click.argument('output_dir', type=click.Path())
@click.option('--augustus_config_path', required=True, help="AUGUSTUS config path")
@click.option('--species', required=True, help="Species name for AUGUSTUS")
def augustus(input_dir, output_dir, augustus_config_path, species):
    """Run AUGUSTUS on fasta files"""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    fasta_files = list(input_dir.glob('*.fa'))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    run_all_augustus(fasta_files, augustus_config_path, species, output_dir)
    merge_augustus_files(output_dir, "augustus.gff3")
    add_exon_to_gene_model("augustus.gff3", "augustus.gff3")

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
    add_exon_to_gene_model("snap.gff3", "snap.gff3")


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
    add_exon_to_gene_model("glimmer.gff3", "glimmer.gff3")


#####################################################
###########         EVM tools         ###############
#####################################################
from wgap.evm_prepare import stringtie_to_gff3, miniprot_to_gene_structure_gff3, miniprot_to_spliced_gff3
from wgap.evm_prepare import prepare_evm_input

@cli.group()
def evm():
    """Subcommand for EVM tools"""




@evm.command('convert')
@click.argument('gff', type=click.Path(exists=True))
@click.option('--source',  default='1', help='input gff source: 1 for stringtie, 2 for miniprot', type=click.Choice(['1', '2']))
@click.option('--miniprot_gff3_format', '-mm', default='1', help='miniprot output GFF3 format: 1 for gene_structure, 2 for spliced_alignment', type=click.Choice(['1', '2']))
@click.option('--output', default='evm.gff3', help='Output file name, default will add evm to the input gff file name')
def format_gff3(gff, source, gff3_format):
    """convert the input gff to gff3 for EVM"""
    # Your format GFF3 logic here
    # 定义映射字典
    source_mapping = {
        '1': 'stringtie',
        '2': 'miniprot', 
    }
    source = source_mapping[source]

    gff3_format_mapping = {
        '1': 'gene_structure',
        '2': 'spliced_alignment'
    }
    gff3_format = gff3_format_mapping[gff3_format]

    if source == 'miniprot':
        out_file = gff.replace('.gff', f'_{gff3_format}_evm.gff3')
    else:
        out_file = gff.replace('.gff', '_evm.gff3')

    if source == "stringtie":
        stringtie_to_gff3(gff, out_file, gff3_format)
    elif source == "miniprot":
        if gff3_format == "gene_structure":
            miniprot_to_gene_structure_gff3(gff, out_file)
        else:
            miniprot_to_spliced_gff3(gff, out_file)
    else:
        log_exception("source not supported")


@evm.command("init")
@click.option('--ab_initio', '-ab', type=str, help='Gene structure from ab initio prediction tools')
@click.option('--transcript','-tr', type=str, help='spliced alignment file from transcriptome alignment')
@click.option('--protein', '-pr', type=str, help='spliced alignment file from transcriptome alignment')
@click.option('--other', '-ot', type=str, help='gene structure file from other')
def init_input_for_evm(ab, other, spliced):
    """Initialize weights for EVM with given parameters."""
    # Split the input strings into lists based on spaces
    ab_initio = ab.split() if ab else []
    transcript = spliced.split() if spliced else []
    protein = other.split() if other else []
    other = other.split() if other else []
    # Your EVM initialization logic here
    prepare_evm_input(ab_initio, transcript, protein, other)


@evm.command('run')
def run_evm():
    """Run EVM"""
    # Your EVM running logic here
    pass






#####################################################
########### post-processing tools     ###############
#####################################################

@cli.group()
def post():
    """"
    post process tools: add_utr, update, rename
    """
    

# add utr with new click
from wgap.utr import add_utr_to_gene_model
@post.command('add_utr')
@click.argument('gff3', type=click.Path(exists=True))
#gtf_files: a list of gtf files
@click.argument('gtf_files', nargs=-1, type=click.Path(exists=True))
@click.option('--output', default='output.gff3', help='Output file name')
def add_utr(gff3, gtf_files, output):
    """Add UTR to GFF3 file"""
    # Your add UTR logic here
    update_genes = add_utr_to_gene_model(gff3, gtf_files)
    # save to gff3
    with open(output, 'w') as f:
        for gene in update_genes.values():
            f.write(gene.to_gff3())
            f.write('\n')


@post.command(
    'rename',
    context_settings = {"ignore_unknown_options" : True},
    short_help = "rename the maker.gff base on given nomenclature"
)
@click.argument(
    'oldgff',
    type = click.Path(dir_okay=True, writable=True, resolve_path=True),
)
@click.argument(
    'newgff',
    type = click.Path(dir_okay=True, writable=True, resolve_path=True),
)
@click.option('-j',
    '--justify',
    help = 'The unique integer portion of the ID will be right justified  with `0`s to this length (default = 5)',
    default = 5,
    type = int
)
@click.option('-p',
    '--prefix',
    type = str,
    help = 'prefix of gene name',
    required=True
)
def run_rename(oldgff, newgff, prefix, justify):
    logging.info("renaming the gff")
    ##############################################
    #  rename the gff
    ##############################################

    # orig_models = maker_update.gff_reader(oldgff)
    # justify = justify
    # out_models = maker_update.update_gene_id(orig_models, prefix, justify, "wgap")
    # maker_update.gff_writer(out_models, newgff)

# # update the protein
@post.command(
    'update',
    context_settings = {"ignore_unknown_options" : True},
    short_help = "update the maker.gff base on the output of apollo"
)
@click.argument(
    'oldgff',
    type = click.Path(dir_okay=True, writable=True, resolve_path=True),
)
@click.argument(
    'newgff',
    type = click.Path(dir_okay=True, writable=True, resolve_path=True),
)
@click.option('-o',
    '--output',
    type = str,
    help = 'output gff',
    required=True
)
def update_gff(oldgff, newgff, outgff):
    logging.info("update the gff")
    # orig_models = maker_update.gff_reader(oldgff)
    # new_models = maker_update.gff_reader(newgff)
    # out_models = maker_update.update_gene_model(orig_models, new_models)
    # maker_update.gff_writer(out_models, outgff)


#####################################################
########### post-processing tools     ###############
#####################################################
from wgap.tools import add_exon_to_gene_model


@cli.group()
def tools():
    """"
    useful tools for wgap
    """

@tools.command(
    'install_ext'
)
def install_ext():
    """
    install the external tools for wgap
    """

# download command for homology evidence
@tools.command(
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

@tools.command(
    'fix'
)
@click.argument(
    "input_file",
    type = click.Path(exists=True)
)
@click.argument(
    "output_file",
    type = click.Path()
)
def add_exon(input_file, output_file):
    """
    add the exon to the gene model
    """
    add_exon_to_gene_model(input_file, output_file)



if __name__ == '__main__':
    cli()
