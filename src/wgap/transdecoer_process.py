from math import inf
import re
import sys
import subprocess

            
def analyze_and_filter_predictions(gene_dict, isoform_dict):
    pattern = r'Name="ORF type:([^,]+) \((\+|\-)\),score=([+-]?\d+\.\d+)"'
    for gene,isoforms in gene_dict.items():
        bad_isoforms = []
        for isofrom in isoforms:
            lines  =  isoform_dict[isofrom]
            predictions = [ line for line in lines if line[2] == 'gene']

            prediction_stat = [re.findall(pattern, line[8])[0] for line in predictions]
            # get the highest score index
            max_score = -inf # very small number
            max_index = -1
            for i, (orf_type, strand, score) in enumerate(prediction_stat):
                if float(score) > max_score:
                    max_score = float(score)
                    max_index = i
            
            # if  the gene with highest score is not compete or strand is not '+' , print
            if max_index == -1:
                print(gene, isofrom, 'no gene prediction')
        
            if prediction_stat[max_index][0] != 'complete' or prediction_stat[max_index][1] != '+':
                #print(gene, isofrom, prediction_stat[max_index][0], prediction_stat[max_index][1])
                bad_isoforms.append(isofrom)
            else:
                best_prediction = predictions[max_index]
                gene_id = best_prediction[8].split(';')[0].split('=')[1]
                # GENE.MSTRG.1.1~~MSTRG.1.1.p1
                prediction_id = gene_id.split('~~')[1]
                newlines = [ best_prediction ]
                for line in lines:
                    if line[2] != 'gene' and prediction_id in line[8]:
                        newlines.append(line)


                #print(gene, isofrom, prediction_stat[max_index][0], prediction_stat[max_index][1])
                
        if len(bad_isoforms) > 0:
            #print(gene, bad_isoforms)
            for isoform in bad_isoforms:
                isoforms.remove(isoform)
        else:
            isoform_dict[isofrom] = newlines
            
    # delete empty gene
    empty_genes = []
    for gene, isoforms in gene_dict.items():
        if len(isoforms) == 0:
            empty_genes.append(gene)

    for gene in empty_genes:
        del gene_dict[gene]
    
    return gene_dict, isoform_dict


def read_gff_file(gff_file):

    gene_dict = {}
    isoform_dict = {}

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip()
            if line == '':
                continue
            line = line.split('\t')

            isoform_id = line[0]
            gene_id = '.'.join(isoform_id.split('.')[0:2])

            if gene_id not in gene_dict:
                gene_dict[gene_id] = set()
            gene_dict[gene_id].add(isoform_id)

            if isoform_id not in isoform_dict:
                isoform_dict[isoform_id] = []

            isoform_dict[isoform_id].append(line)
    
    return gene_dict, isoform_dict

# write the result in transdecoder.gff3
def write_gff3(gene_dict, isoform_dict, out_file):
    with open(out_file, 'w') as f:
        for gene, isoforms in gene_dict.items():
            for isoform in isoforms:
                lines = isoform_dict[isoform]
                for line in lines:
                    f.write('\t'.join(line) + '\n')
                f.write('\n')


def main(gtf_file, genome_file, prefix):

    cnda_file = f"{prefix}.fasta"
    cmd = f"gtf_genome_to_cdna_fasta {gtf_file} {genome_file} > {cnda_file}"
    subprocess.run(cmd, shell=True, check=True)

    # predict the orf
    cmd = f"TransDecoder.LongOrfs -t {cnda_file}"
    subprocess.run(cmd, shell=True, check=True)

    # predict the orf
    cmd = f"TransDecoder.Predict -t {cnda_file}"
    subprocess.run(cmd, shell=True, check=True)

    gff_file = f"{cnda_file}.transdecoder.gff3"
    filtered_file = f"{cnda_file}.transdecoder.filtered.gff3"

    # read the gff file
    gene_dict, isoform_dict = read_gff_file(gff_file)
    gene_dict, isoform_dict = analyze_and_filter_predictions(gene_dict, isoform_dict)
    write_gff3(gene_dict, isoform_dict, filtered_file)


    cmd = f"gtf_to_alignment_gff3.pl {gtf_file} > {prefix}.genome.gtff3"
    subprocess.run(cmd, shell=True, check=True)

    cmd = f"cdna_alignment_orf_to_genome_orf.pl {filtered_file} {prefix}.genome.gtff3 {cnda_file} > {prefix}.cdna.genome.gff3"
    subprocess.run(cmd, shell=True, check=True)


if __name__ == '__main__':
    gtf_file = sys.argv[1]
    genome_file = sys.argv[2]
    prefix= sys.argv[3]
    main(gtf_file, genome_file, prefix)

