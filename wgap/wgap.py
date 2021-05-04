import os
import sys
import click
import logging
import subprocess

from shutil import copyfile
from snakemake import load_configfile

from wgap.scripts import utils

from wgap import __version__

#from wagp.conf import run_init


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

# define command line
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """WGAP - a whole genome annotation pipeline
    """

#cli.add_command(run_init)


# workflow command
@cli.command(
    "run",
    context_settings = {"ignore_unknown_options" : True},
    short_help = "run wgap main workflow"
)
@click.argument(
    "workflow",
    default="all",
    type=click.Choice(["None", "all", "test"]),
)
@click.option("-w",
    "--working-dir",
    type = click.Path(dir_okay = True, writable=True, resolve_path=True),
    help = "project running directory",
    default = "."
)
@click.option("-c",
    "--config-file",
    type = click.Path(dir_okay = True, writable=True, resolve_path=True),
    help = "config file, generated with wgap init",
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    help="use at most this many jobs in parallel (see cluster submission for mor details).",
)
@click.option(
    "--profile",
    default=None,
    help="snakemake profile e.g. for cluster execution.",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(workflow, working_dir, config_file, jobs, profile, dryrun, snakemake_args):

    # check file availability
    if config_file is None:
        config_file = os.path.join(working_dir, 'config.yaml')
    
    if not os.path.exists(config_file):
        logging.critical(f'config-file not found: {config_file}.')
        sys.exit(1)
    
    conf = load_configfile(config_file)
    sample_file = conf['sample']

    if not os.path.exists(sample_file):
        logging.critical('{} is not existed'.format(sample_file))
        sys.exit(1)
    
    # TO DO: check file avaiable
    
    cmd = (
        "snakemake --snakefile {snakefile} --directory {working_dir}"
        " {jobs} --rerun-incomplete "
        " --configfile '{config_file}' --nolock "
        " {profile} {dryrun}"
        " {target_rule} "
        " {args} "
    ).format(
        snakefile = utils.get_snakefile(),
        working_dir = working_dir,
        jobs = "--jobs {}".format(jobs if jobs is not None else "4"),
        config_file = config_file,
        profile = "" if (profile is None) else "--profile {}".foramt(profile),
        dryrun="--dryrun" if dryrun else "",
        args = " ".join(snakemake_args),
        target_rule = workflow if workflow != "None" else ""
    )
    logging.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logging.critical(e)
        exit(1)

@cli.command(
    'init',
    context_settings = {"ignore_unknown_options" : True},
    short_help = "copy the config.yaml and samples.csv to working directoty"
)
@click.argument(
    'workdir',
    default="."
)
def init_workdir(workdir):
    config_file = utils.get_configfile()
    sample_file = utils.get_samplefile()

    copyfile(config_file, os.path.join(workdir, "config.yaml"))
    copyfile(sample_file, os.path.join(workdir, "samples.csv"))
    

# download command for homology evidence
@cli.command(
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
    dat_file = utils.download_uniprot(specie, dataset)
    if not dat_file:
        sys.exit(1)
    utils.convet_dat_to_fasta(dat_file, fasta)
    logging.info("Finished: %s" % fasta )



if __name__ == "__main__":
    cli()