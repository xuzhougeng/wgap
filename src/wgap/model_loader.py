
from .tools import extract_chromosome_coordinates
from typing import Dict, List

from .gene import Gene, Transcript, Exon, CDS
from .seq_utils import SeqRecordDict

def get_attributes(attributes: str) -> Dict[str, str]:
    """
    Parses the attributes column of a GTF file and returns a dictionary of attributes.

    Args:
        attributes (str): The attributes column of a GTF file.

    Returns:
        Dict[str, str]: A dictionary of attributes.
    """
    attributes_dict = {}
    for attribute in attributes.split(";"):
        if attribute.strip() == "":
            continue

        key, value = attribute.strip().split(" ")
        attributes_dict[key] = value.replace('"', '')
    return attributes_dict

def stringtie_loader(stringtie_gtf: str, fasta_dict: SeqRecordDict = None) -> Dict[str, Gene]:
    
    gene_dict : Dict[str, Gene] = {} 
    exon_dict : Dict[str, Exon] = {}
    
    for line in open(stringtie_gtf):
        if line.startswith("#"):
            continue
        line = line.strip()

        # Split the line into columns
        fields = line.split("\t")
        if len(fields) != 9:
            continue

        chrom = fields[0]
        feature_type = fields[2]

        # gff is 1-based, convert to 0-based
        start, end = int(fields[3]) -1 , int(fields[4])
        strand = fields[6]
        attributes = fields[8]
        attributes_dict = get_attributes(attributes)
        
        if feature_type == "transcript":
            
            gene_id = attributes_dict["gene_id"]
            transcript_id = attributes_dict["transcript_id"]

            if gene_id not in gene_dict:
                gene_dict[gene_id] = Gene(gene_id)
            
            transcript = Transcript(transcript_id, chrom, start, end, strand)
            gene_dict[gene_id].add_transcript( transcript )

        elif feature_type == "exon":

            sequence = fasta_dict[chrom][start:end]
            exon_id = f"{chrom}_{start}_{end}"

            if exon_id not in exon_dict:
                exon_dict[exon_id] = Exon(exon_id, chrom, start, end, strand, sequence)

            # add exon to transcript
            gene_id = attributes_dict["gene_id"]
            transcript_id = attributes_dict["transcript_id"]
            
            transcript = gene_dict[gene_id].get_transcript(transcript_id)
            transcript.add_exon(exon_dict[exon_id])
        else:
            continue

    
    return gene_dict

def augustus_loader(augustus_gtf: str, fasta_dict: SeqRecordDict = None) -> Dict[str, Gene]:

    gene_dict: Dict[str, Gene] = {}
    cds_dict: Dict[str, CDS] = {}

    file_handle = open(augustus_gtf)
    for line in file_handle:
        if line.startswith("#"):
            continue
        line = line.strip()
        
        parts = line.split("\t")
        if len(parts) != 9:
            continue  # Ignore malformed lines
       
        chrom, _, feature_type, start, end, _, strand, phase, attributes = parts
        
        if len(chrom.split("_")) > 1:
            chrom, sub_start, _ = extract_chromosome_coordinates(chrom)
            # gff is 1-based, convert to 0-based
            start, end = int(start) - 1 + int(sub_start), int(end) + int(sub_start)
        else:
            # gff is 1-based, convert to 0-based
            start, end = int(start) - 1, int(end)

        if feature_type == "gene":
            gene_id = attributes
            new_gene_id = f"{chrom}_{start}_{end}_{gene_id}"
            gene = Gene(gene_id=new_gene_id, chrom=chrom, start=start, end=end, strand=strand)
            # add gene to gene_dict
            if gene_id not in gene_dict:
                gene_dict[gene_id] = gene

        elif feature_type == "transcript":
            tx_id = f"{gene.id}_{attributes}"
            gene_id = attributes.split(".")[0]
            transcript = Transcript(transcript_id=tx_id, chrom=chrom, start=start, end=end, strand=strand)
            gene = gene_dict[gene_id]
            gene.add_transcript(transcript)

        elif feature_type == "CDS":
            cds_id = f"CDS_{chrom}_{start}_{end}"
            #sequence = fasta_dict[chrom][start:end]  # Get the exon sequence from the genome
            cds = CDS(cds_id=cds_id, chrom=chrom, start=start, end=end, strand=strand, phase=int(phase))
            if fasta_dict is not None:
                cds.sequence = fasta_dict[chrom][start:end]
            if cds_id not in cds_dict:
                cds_dict[cds_id] = cds

            tx_id = get_attributes(attributes)["transcript_id"]
            gene_id = tx_id.split(".")[0]
            gene = gene_dict[gene_id]

            tx_key =f"{gene.id}_{tx_id}"
            transcript = gene.get_transcript(tx_key)
            
            if transcript is not None:
                transcript.add_cds(cds)

    file_handle.close()

    return gene_dict

def gff3_loader(gff3_file: str, fasta_dict: SeqRecordDict = None) -> Dict[str, Gene]:
    
    gene_dict : Dict[str, Gene] = {} 
    exon_dict : Dict[str, Exon] = {}
    cds_dict : Dict[str, CDS] = {}
    
    file_handle = open(gff3_file)
    for line in file_handle:
        if line.startswith("#"):
            continue
        line = line.strip()

        #parts = re.split(r'\s+', line.strip())
        parts = line.split("\t")
        if len(parts) != 9:
            continue  # Ignore malformed lines

        chrom, _, feature_type, start, end, _, strand, phase, attributes = parts
        # gff is 1-based, convert to 0-based
        start, end = int(start) - 1, int(end)

        attributes_dict = dict(item.split("=") for item in attributes.split(";") if "=" in item)

        # check fasta_dict is a dict
        if  type(fasta_dict) is  dict and chrom not in fasta_dict:
            continue

        if feature_type == "gene":
            gene_id = attributes_dict["ID"]
            gene = Gene(gene_id=gene_id, chrom=chrom, start=start, end=end, strand=strand)
            if gene_id not in gene_dict:
                gene_dict[gene_id] = gene

        elif feature_type == "mRNA":
            transcript_id = attributes_dict["ID"]
            gene_id = attributes_dict["Parent"]
            transcript = Transcript(transcript_id=transcript_id, chrom=chrom, start=start, end=end, strand=strand)
            gene = gene_dict[gene_id]
            gene.add_transcript(transcript)

        elif feature_type == "exon":
            exon_id = attributes_dict.get("ID", f"exon_{chrom}_{start}_{end}")
            transcript_id = attributes_dict["Parent"]

            sequence = None
            if fasta_dict is not None:
                sequence = fasta_dict[chrom][start:end]  # Get the exon sequence from the genome

            exon = Exon(exon_id=exon_id, chrom=chrom, start=start, end=end, strand=strand, sequence=sequence)
            if exon.id not in exon_dict:
                exon_dict[exon.id] = exon
            transcript = gene_dict[gene_id].get_transcript(transcript_id)
            if transcript is not None:
                transcript.add_exon(exon)
        elif feature_type == "CDS":
            cds_id = attributes_dict.get("ID", f"cds_{chrom}_{start}_{end}")
            transcript_id = attributes_dict["Parent"]
            
            sequence = None
            if fasta_dict is not None:
                sequence = fasta_dict[chrom][start:end]  # Get the exon sequence from the genome
            
            cds = CDS(cds_id=cds_id, chrom=chrom, start=start, end=end, strand=strand, sequence=sequence, phase=phase)
            if cds.id not in exon_dict:
                cds_dict[cds.id] = cds
            transcript = gene_dict[gene_id].get_transcript(transcript_id)
            if transcript is not None:
                transcript.add_cds(cds)
    file_handle.close()
    
    return gene_dict

def zff_loader(zff_file_path: str, fasta_dict: SeqRecordDict = None) -> Dict[str, Gene]:
    """
    ZFF format long
    Label, Begin, End,  Strand, Score, 5’-overhang[phase in GFF], 3’-overhang, Frame, Group
    5’- and 3’-overhang are the number of bp of an incomplete codon at each end of an exon
    """
    gene_dict: Dict[str, Gene] = {}
    cds_dict: Dict[str, List[CDS]] = {}

    with open(zff_file_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                continue  # Skip header lines
            line = line.strip()

            # Process feature lines
            parts = line.split('\t')
            if len(parts) != 9:
                continue  # Ensure line has enough parts for parsing
            
            feature_type, start, end, strand, score, phase, _, _, attributes = parts
            start, end = int(start), int(end)
            
            # attributes format is : chrom-snap.id
            coord, gene_id = attributes.split("-")
            
            if extract_chromosome_coordinates(coord) is not None:
                chrom, rel_start, _ = extract_chromosome_coordinates(coord)
                start = start + rel_start
                end = end + rel_start
            else:
                chrom = coord
            
            cds_id = f"CDS_{chrom}_{start}_{end}"
            cds = CDS(cds_id=cds_id, chrom=chrom, start=start, end=end, strand=strand, phase=int(phase))
            if fasta_dict is not None:
                cds.sequence = fasta_dict[chrom][start:end]

            if gene_id not in cds_dict:
                cds_dict[gene_id] = [ cds ]
            else:
                cds_dict[gene_id].append(cds)
        
    for gene_id, cds_list in cds_dict.items():
        new_gene_id = f"{chrom}_{start}_{end}_{gene_id}"
        tx_id = f"{new_gene_id}.t1"
        # sort cds_list

        cds_list = sorted(cds_list, key=lambda cds: cds.start)
        gene_start = cds_list[0].start
        gene_end = cds_list[-1].end
        gene_strand = cds_list[0].strand
        gene = Gene(gene_id=new_gene_id, chrom=chrom, start=gene_start, end=gene_end, strand=gene_strand)

        tx = Transcript(transcript_id=tx_id, chrom=chrom, start=gene_start, end=gene_end, strand=gene_strand)
        
        for cds in cds_list:
            tx.add_cds(cds)
        gene.add_transcript(tx)

        gene_dict[gene_id] = gene
        
    return gene_dict

def glimmer_loader(glimmer_file_path: str, fasta_dict: SeqRecordDict = None) -> Dict[str, Gene]:
    gene_dict: Dict[str, Gene] = {}
    cds_dict: Dict[str, List[CDS]] = {}
    
    # 读取所有行
    glimmer_result = open(glimmer_file_path).readlines()
    # 删除空行
    glimmer_result = [x.strip() for x in glimmer_result if x.strip() != ""]

    # 解析序列名
    chrom = ""
    chrom_start = 0
    sequence_name = glimmer_result[1][15:]
    coord = extract_chromosome_coordinates(sequence_name)
    if coord is not None:
        chrom = coord[0]
        chrom_start = int(coord[1])
    else:
        chrom = sequence_name

    # 解析结果
    for line in glimmer_result[6:]:
        parts = line.split()
        if len(parts) != 7:
            continue  # Skip malformed lines
        
        gene_id, exon_id, strand, exon_type, start, end, exon_length = parts
        start, end = int(start) + chrom_start, int(end) + chrom_start  # Convert start and end to integers
        
        cds_id = f"CDS_{chrom}_{start}_{end}"
        cds = CDS(cds_id=cds_id, chrom=chrom, start=start, end=end, strand=strand)
        if fasta_dict is not None:
            cds.sequence = fasta_dict[chrom][start:end]

        gene_id =   "g" + gene_id

        if gene_id not in cds_dict:
            cds_dict[gene_id] = [ cds ]
        else:
            cds_dict[gene_id].append(cds)
        
    for gene_id, cds_list in cds_dict.items():
        new_gene_id = f"{chrom}_{start}_{end}_{gene_id}"
        tx_id = f"{new_gene_id}.t1"
        # sort cds_list

        cds_list = sorted(cds_list, key=lambda cds: cds.start)
        gene_start = cds_list[0].start
        gene_end = cds_list[-1].end
        gene_strand = cds_list[0].strand
        gene = Gene(gene_id=new_gene_id, chrom=chrom, start=gene_start, end=gene_end, strand=gene_strand)

        tx = Transcript(transcript_id=tx_id, chrom=chrom, start=gene_start, end=gene_end, strand=gene_strand)
        
        for cds in cds_list:
            tx.add_cds(cds)
        gene.add_transcript(tx)

        gene_dict[gene_id] = gene    
    
    return gene_dict
