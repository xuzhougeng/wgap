
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from .gene import gene_to_gff3
from .model_loader import augustus_loader, zff_loader, glimmer_loader

MAX_WORKERS = 10

def run_augustus(fasta_file, AUGUSTUS_CONFIG_PATH, species, 
                 output_dir="augustus_dir") -> None:
    base_name = fasta_file.stem
    output_file = output_dir / f"{base_name}.gff3"
    status_file = output_dir / f"{base_name}.augustus.done"

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
    output_file = output_dir / f"{base_name}.zff"
    status_file = output_dir / f"{base_name}.snap.done"

    # Skip if task already completed
    if status_file.exists():
        print(f"Skipping {fasta_file}, already processed.")
        return
    
    with open(output_file, 'w') as f_out:
        cmd = ["snap", f"{snap_hmm_file}", str(fasta_file)]
        result = subprocess.run(cmd, stdout=f_out, stderr=subprocess.PIPE)

    if result.returncode == 0:
        # Write a status file on successful completion
        status_file.touch()
    else:
        raise RuntimeError(f"SNAP failed on {fasta_file}: {result.stderr}")

def run_glimmer(fasta_file, glimmerhmm_model_dir, output_dir="glimmer_out") -> None:
    base_name = fasta_file.stem

    output_file = output_dir / f"{base_name}.glimmer"
    status_file = output_dir / f"{base_name}.glimmer.done"

    # Skip if task already completed
    if status_file.exists():
        print(f"Skipping {fasta_file}, already processed.")
        return

    
    with open(output_file, 'w') as f_out:
        cmd = ["glimmerhmm", str(fasta_file), glimmerhmm_model_dir]
        result = subprocess.run(cmd, stdout=f_out, stderr=subprocess.PIPE)

    if result.returncode == 0:
        # Write a status file on successful completion
        status_file.touch()
    else:
        raise RuntimeError(f"GLIMMER failed on {fasta_file}: {result.stderr}")

def parallel_process(function, fasta_files, *args):
    """
    并行执行指定的函数 `function`，对 `fasta_files` 列表中的每个文件应用该函数。
    `*args` 是传递给 `function` 的其他参数。
    """
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(function, fasta_file, *args) for fasta_file in fasta_files]
        for future in futures:
            future.result()  # 等待所有任务完成

def run_all_augustus(fasta_files, AUGUSTUS_CONFIG_PATH, species, output_dir):
    """
    并行执行 AUGUSTUS。
    """
    parallel_process(run_augustus, fasta_files, AUGUSTUS_CONFIG_PATH, species, output_dir)

def run_all_snap(fasta_files, snap_hmm_file, output_dir):
    """
    并行执行 SNAP。
    """
    parallel_process(run_snap, fasta_files, snap_hmm_file, output_dir)

def run_all_glimmer(fasta_files, glimmerhmm_model_dir, output_dir):
    """
    并行执行 GLIMMER。
    """
    parallel_process(run_glimmer, fasta_files, glimmerhmm_model_dir, output_dir)

def merge_augustus_files(input_dir, output_file) -> None:
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

def merge_snap_files(input_dir, output_file) -> None:
    with open(output_file, 'w') as outfile:
        # Iterate over each zff file in the input directory
        for zff_file in os.listdir(input_dir):
            if zff_file.endswith(".zff"):
                full_path = os.path.join(input_dir, zff_file)
                
                # Load the zff file using the provided loader
                genes = zff_loader(full_path)
                if len(genes) == 0:
                    continue

                # Convert and write each gene to the output file
                for gene_id, gene in genes.items():
                    gff3_output = gene_to_gff3(gene)
                    outfile.write(gff3_output + "\n")

def merge_glimmer_files(input_dir, output_file) -> None:
    with open(output_file, 'w') as outfile:
        # Iterate over each glimmer file in the input directory
        for glimmer_file in os.listdir(input_dir):
            if glimmer_file.endswith(".glimmer"):
                full_path = os.path.join(input_dir, glimmer_file)
                
                # Load the glimmer file using the provided loader
                genes = glimmer_loader(full_path)
                if len(genes) == 0:
                    continue

                # Convert and write each gene to the output file
                for gene_id, gene in genes.items():
                    gff3_output = gene_to_gff3(gene)
                    outfile.write(gff3_output + "\n")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Gene prediction tool wrapper")
    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # AUGUSTUS 子命令
    parser_augustus = subparsers.add_parser('augustus', help='Run AUGUSTUS on fasta files')
    parser_augustus.add_argument("input_dir", type=Path, help="Directory containing fasta files")
    parser_augustus.add_argument("output_dir", type=Path, help="Output directory for AUGUSTUS gff3 files")
    parser_augustus.add_argument("--AUGUSTUS_CONFIG_PATH", required=True, help="AUGUSTUS config path")
    parser_augustus.add_argument("--species", required=True, help="Species name for AUGUSTUS")

    # SNAP 子命令
    parser_snap = subparsers.add_parser('snap', help='Run SNAP on fasta files')
    parser_snap.add_argument("input_dir", type=Path, help="Directory containing fasta files")
    parser_snap.add_argument("output_dir", type=Path, help="Output directory for SNAP output files")
    parser_snap.add_argument("--snap_hmm_file", required=True, help="HMM file for SNAP")

    # GLIMMER 子命令
    parser_glimmer = subparsers.add_parser('glimmer', help='Run GlimmerHMM on fasta files')
    parser_glimmer.add_argument("input_dir", type=Path, help="Directory containing fasta files")
    parser_glimmer.add_argument("output_dir", type=Path, help="Output directory for GlimmerHMM output files")
    parser_glimmer.add_argument("--glimmerhmm_model_dir", required=True, help="Model directory for GLIMMER")

    args = parser.parse_args()

    fasta_files = list(args.input_dir.glob('*.fa'))  # 获取输入目录下的所有 .fa 文件
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    if args.command == 'augustus':
        run_all_augustus(fasta_files, args.AUGUSTUS_CONFIG_PATH, args.species, args.output_dir)
        merge_augustus_files(args.output_dir, "augustus.gff3")
    elif args.command == 'snap':
        run_all_snap(fasta_files, args.snap_hmm_file, args.output_dir)
        merge_snap_files(args.output_dir, "snap.gff3")
    elif args.command == 'glimmer':
        run_all_glimmer(fasta_files, args.glimmerhmm_model_dir, args.output_dir)
        merge_glimmer_files(args.output_dir, "glimmer.gff3")
    else:
        parser.print_help()

if __name__ == "__main__":
    main()