
import re
import logging

logging.basicConfig(level=logging.INFO)


def stringtie_to_gff3(stringtie_gtf, gff3_out, source="StringTe"):
    """
    Convert stringtie gtf to EVM spliced alignment gff3 formatting
    """
    with open(stringtie_gtf, "r") as f, open(gff3_out, "w") as evmfile:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            chrom, _, feature, start, end, score, strand, frame, attributes = parts
            if feature == "transcript":
                continue
            # get transcript id in attributes
            match = re.search(r'transcript_id "([^"]+)"', attributes)
            tx_id = match.group(1)

            feature = "EST_match"

            out_line = f"{chrom}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={tx_id}\n"
            evmfile.write(out_line)

def miniprot_to_spliced_gff3(miniprot_gtf, gff3_out, source="miniprot_protAln"):
    """
    Reimplement the function from  miniprot_GFF_2_EVM_alignment_GFF3.py which contributed by: Kresimir.Krizanovic@fer.hr
    """
    with open(miniprot_gtf, "r") as csvfile, open(gff3_out, "w") as evmfile:
        # Read the first line
        line = csvfile.readline()
        i = 1

        # Structures for storing a data for one gene/transcript
        # miniprot currently contains lines for mrna, CDS and stop_codon records
        mrna_line = []
        cds_lines = []  # A list of CDS lines
        sc_line = []

        start = True
        while line:
            # This is just in case
            if line == "":
                break

            if line.startswith(
                "##PAF"
            ):  # transform all previously collected data and write it into the output file
                if (
                    not start and len(mrna_line) > 0
                ):  # Just not at the start, because data has not been collected
                    # no gene line for spliced alignments
                    mrna_params = mrna_line[8]
                    pos = mrna_params.find(";")
                    mrna_id = mrna_params[3:pos]

                    # CDS records do not have an ID, for spliced alignments, all CDS records have the same ID
                    # Using mRNa ID
                    cnt = 0
                    for cds_line in cds_lines:
                        cnt += 1
                        # no exon lines for spliced alignments
                        cds_line[1] = source
                        cds_line[2] = "nucleotide_to_protein_match"
                        # modifying CDS line params to set ID and setting score fiels as identity percentage
                        cds_elems = cds_line[8].split(";")
                        mrna_elems = mrna_line[8].split(";")
                        identity = float(cds_elems[2][9:])  # Extracting identity parameter
                        cds_line[5] = "{0:.2f}".format(
                            identity * 100
                        )  # Setting identitty percentage as score
                        cds_elems[0] = mrna_elems[0]  # Using MRNA ID for all CDS alignments
                        cds_line[8] = ";".join(cds_elems)

                        evmfile.write("\t".join(cds_line))

                    mrna_line = []
                    cds_lines = []  # A list of CDS lines
                    sc_line = []

            if line.startswith("#"):  # Comment line, just write it in the output file
                # evmfile.write(line)
                pass
            else:
                elements = line.split("\t")
                if elements[2].upper() == "MRNA":
                    start = False  # We have collected some data - not start any more
                    mrna_line = elements
                if elements[2].upper() == "CDS":
                    cds_lines.append(elements)

                # Stop_codon usually marks the end of the transcript, but not always

                if elements[2].upper() == "STOP_CODON":
                    sc_line = elements

            line = csvfile.readline()
            i += 1

        # Write the last collected set of data
        mrna_params = mrna_line[8]
        pos = mrna_params.find(";")
        mrna_id = mrna_params[3:pos]

        # Using mRNa ID
        cnt = 0
        for cds_line in cds_lines:
            cnt += 1
            # no exon lines for spliced alignments
            cds_line[1] =  source
            cds_line[2] = "nucleotide_to_protein_match"
            # modifying CDS line params to set ID and setting score fiels as identity percentage
            cds_elems = cds_line[8].split(";")
            mrna_elems = mrna_line[8].split(";")
            identity = float(cds_elems[2][9:])  # Extracting identity parameter
            cds_line[5] = "{0:.2f}".format(
                identity * 100
            )  # Setting identitty percentage as score
            cds_elems[0] = mrna_elems[0]  # Using MRNA ID for all CDS alignments
            cds_line[8] = ";".join(cds_elems)

            evmfile.write("\t".join(cds_line))

    logging.info(f"Done! Read {i} lines")

def miniprot_to_gene_structure_gff3(miniprot_gtf, gff3_out, source="MiniProt"):
    """
    Reimplement the function from  miniprot_GFF_2_EVM_alignment_GFF3.py which contributed by: Kresimir.Krizanovic@fer.hr
    """
    with open(miniprot_gtf, "r") as csvfile, open(gff3_out, "w") as evmfile:
        # Read the first line, should be a comment
        line = csvfile.readline()
        i = 1

        # Structures for storing a data for one gene/transcript
        # miniprot currently contains lines for mrna, CDS and stop_codon records
        gene_line = []
        mrna_line = []
        cds_lines = []  # A list of CDS lines
        sc_line = []

        start = True
        while line:
            # This is just in case
            if line == "":
                break

            if line.startswith(
                "##PAF"
            ):  # transform all previously collected data and write it into the output file
                if (
                    not start and len(mrna_line) > 0
                ):  # Just not at the start, because data has not been collected
                    gene_line = mrna_line[:]  # copy the mrna data
                    gene_line[2] = "gene"
                    mrna_params = mrna_line[8]
                    pos = mrna_params.find(";")
                    mrna_id = mrna_params[3:pos]
                    gene_id = "G_" + mrna_id
                    gene_line[8] = "ID={0};Name={1} model {2}\n".format(
                        gene_id, "miniprot", mrna_id
                    )
                    mrna_line[8] = (
                        mrna_params[: pos + 1] + "Parent=" + gene_id + mrna_params[pos:]
                    )

                    evmfile.write("\t".join(gene_line))
                    evmfile.write("\t".join(mrna_line))

                    # CDS records do not have an ID, have to make one
                    cnt = 0
                    for cds_line in cds_lines:
                        cnt += 1
                        exon_line = cds_line[:]  # copy all elements
                        exon_line[2] = "exon"
                        exon_line[7] = "."  # set phase field for exons to '.'
                        exon_line[8] = (
                            "ID={0}.{1};".format(mrna_id, cnt) + exon_line[8]
                        )  # set exon id
                        cds_line[8] = (
                            "ID=CDS_{0}.{1};".format(mrna_id, cnt) + cds_line[8]
                        )  # set cds id

                        evmfile.write("\t".join(exon_line))
                        evmfile.write("\t".join(cds_line))

                    if len(sc_line) > 0:
                        # Skip stop codon line
                        # evmfile.write('##' + '\t'.join(sc_line))
                        pass

                    gene_line = []
                    mrna_line = []
                    cds_lines = []  # A list of CDS lines
                    sc_line = []

            if line.startswith("#"):  # Comment line, just skip
                # evmfile.write(line)
                pass
            else:
                elements = line.split("\t")
                elements[1] = source
                if elements[2].upper() == "MRNA":
                    start = False  # We have collected some data - not start any more
                    mrna_line = elements
                if elements[2].upper() == "CDS":
                    cds_lines.append(elements)

                # Stop_codon usually marks the end of the transcript, but not always

                if elements[2].upper() == "STOP_CODON":
                    sc_line = elements

            line = csvfile.readline()
            i += 1

        # Write the last collected set of data
        gene_line = mrna_line[:]  # copy the mrna data
        gene_line[2] = "gene"
        mrna_params = mrna_line[8]
        pos = mrna_params.find(";")
        mrna_id = mrna_params[3:pos]
        gene_id = "G_" + mrna_id
        gene_line[8] = "ID={0};Name={1} model {2}\n".format(gene_id, "miniprot", mrna_id)
        mrna_line[8] = mrna_params[: pos + 1] + "Parent=" + gene_id + mrna_params[pos:]

        evmfile.write("\t".join(gene_line))
        evmfile.write("\t".join(mrna_line))

        # CDS records do not have an ID, have to make one
        cnt = 0
        for cds_line in cds_lines:
            cnt += 1
            exon_line = cds_line[:]  # copy all elements
            exon_line[2] = "exon"
            exon_line[7] = "."  # set phase field for exons to '.'
            exon_line[8] = "ID={0}.{1};".format(mrna_id, cnt) + exon_line[8]  # set exon id
            cds_line[8] = "ID=CDS_{0}.{1};".format(mrna_id, cnt) + cds_line[8]  # set cds id

            evmfile.write("\t".join(exon_line))
            evmfile.write("\t".join(cds_line))

        if len(sc_line) > 0:
            # skip stop codon line
            # evmfile.write('##' + '\t'.join(sc_line))
            pass

    logging.info(f"Done! Read {i} lines")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Convert gtf/gff to EVM gf3")
    parser.add_argument("input", help="Input gtf/gff file")
    parser.add_argument("output", help="Output gff3 file")
    parser.add_argument("--source", default="StringTie", help="Source of the input file")
    parser.add_argument("--type", default="stringtie", choices=["stringtie", "miniprot"], help="Type of the input file")
    parser.add_argument("--output_type", default="spliced", choices=["spliced", "gene_structure"], help="Type of the output file")
    args = parser.parse_args()

    if args.type == "stringtie":
        stringtie_to_gff3(args.input, args.output, args.source)
    elif args.type == "miniprot":
        if args.output_type == "spliced":
            miniprot_to_spliced_gff3(args.input, args.output, args.source)
        elif args.output_type == "gene_structure":
            miniprot_to_gene_structure_gff3(args.input, args.output, args.source)
    else:
        logging.error("Unknown input type")
        return 1
    


if __name__ == "__main__":
    main()