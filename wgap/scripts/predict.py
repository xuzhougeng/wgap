
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from wgap.scripts.gene import augustus_loader, gene_to_gff3

MAX_WORKERS = 10

def run_augustus(fasta_file, AUGUSTUS_CONFIG_PATH, species, 
                 output_dir="augustus_dir") -> None:
    base_name = fasta_file.stem
    output_file = output_dir / f"{base_name}.gff3"
    status_file = output_dir / f"{base_name}.status"

    # Skip if task already completed
    if status_file.exists():
        print(f"Skipping {fasta_file}, already processed.")
        return

    cmd = ["augustus", f"--AUGUSTUS_CONFIG_PATH={AUGUSTUS_CONFIG_PATH}", f"--species={species}", str(fasta_file), f"--outfile={output_file}"]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode == 0:
        # Write a status file on successful completion
        status_file.touch()
    else:
        raise RuntimeError(f"AUGUSTUS failed on {fasta_file}: {result.stderr}")


def run_snap(fasta_file, snap_hmm_file, output_dir="snap_dir") -> None:
    base_name = fasta_file.stem
    output_file = output_dir / f"{base_name}.gff3"
    status_file = output_dir / f"{base_name}.status"

    # Skip if task already completed
    if status_file.exists():
        print(f"Skipping {fasta_file}, already processed.")
        return

    cmd = ["snap", f"{snap_hmm_file}", str(fasta_file), "-o", f"{output_file}"]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode == 0:
        # Write a status file on successful completion
        status_file.touch()
    else:
        raise RuntimeError(f"SNAP failed on {fasta_file}: {result.stderr}")

def run_glimmer(fasta_file, glimmerhmm_model_dir, output_dir="glimmer_out") -> None:
    base_name = fasta_file.stem

    output_file = output_dir / f"{base_name}.gff3"
    status_file = output_dir / f"{base_name}.status"

    # Skip if task already completed
    if status_file.exists():
        print(f"Skipping {fasta_file}, already processed.")
        return
    # glimmerhmm 01_35295998_35318392.fa ol
    cmd = ["glimmerhmm", str(fasta_file), glimmerhmm_model_dir]

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode == 0:
        # Write a status file on successful completion
        status_file.touch()
    else:
        raise RuntimeError(f"GLIMMER failed on {fasta_file}: {result.stderr}")


def merge_gff3_files(input_dir, output_file) -> None:
    with open(output_file, 'w') as outfile:
        # Iterate over each gff3 file in the input directory
        for gff3_file in os.listdir(input_dir):
            if gff3_file.endswith(".gff3"):
                full_path = os.path.join(input_dir, gff3_file)
                
                # Load the gff3 file using the provided loader
                genes = augustus_loader(full_path)
                if len(genes) == 0:
                    continue

                # Convert and write each gene to the output file
                for gene_id, gene in genes.items():
                    gff3_output = gene_to_gff3(gene)
                    outfile.write(gff3_output + "\n")

def main(input_dir, output_dir, AUGUSTUS_CONFIG_PATH, species): 
    
    # Create a pool of workers and execute the jobs
    with ThreadPoolExecutor(max_workers = MAX_WORKERS) as executor:
        futures = [executor.submit(run_augustus, fasta_file, output_dir, AUGUSTUS_CONFIG_PATH, species) for fasta_file in input_dir.glob("*.fa")]
        for future in futures:
            future.result()  # Wait for all futures to complete

    # Merge the gff3 files
    output_file = str( Path(".") / "augustus.gff3")
    merge_gff3_files(output_dir, output_file)

if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description="Run AUGUSTUS on a set of fasta files")
    parser.add_argument("input_dir", help="Directory containing fasta files")
    parser.add_argument("output_dir", help="Output directory for gff3 files")
    # species
    parser.add_argument("-s", "--species", help="species")
    # config path
    parser.add_argument("-c", "--AUGUSTUS_CONFIG_PATH", help="AUGUSTUS_CONFIG_PATH")

    args = parser.parse_args()

    config_path = args.AUGUSTUS_CONFIG_PATH 

    if config_path is None:
        config_path = os.environ.get("AUGUSTUS_CONFIG_PATH")

    if config_path is None or len(config_path) == 0:
        raise RuntimeError("AUGUSTUS_CONFIG_PATH not provided")
    if args.species is None:
        raise RuntimeError("species not provided")
    
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    species = args.species
    

    # Create a pool of workers and execute the jobs
    main(input_dir, output_dir, config_path, species)

