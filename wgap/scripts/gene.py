from typing import Optional, Dict, List
from wgap.scripts.seq_utils import reverse_complement, translate_cds, SeqRecordDict
from wgap.scripts.tools import extract_chromosome_coordinates

# 基因结构: gene->(mRNA->(exon->cds(?))(+) )(+) 
# 一个基因可以有多个mRNA, 一个mRNA可以有多个exon, 一个外显子最多只有一个cds, 但是如果完全是非翻译的(UTR), 则可能没有cds片段

class Exon:
    """
    Represents an exon in a gene.

    Attributes:
        id (str): The ID of the exon.
        chrom (str): The chromosome where the exon is located.
        start (int): The start position of the exon.
        end (int): The end position of the exon.
        strand (str): The strand of the exon. "-" or "+"
        sequence (str): The sequence of the exon, same with genome.
        size (int): The size of the exon sequence.

    Methods:
        __repr__(): Returns a string representation of the exon.
    """

    def __init__(self, exon_id: Optional[str] = None, 
                 chrom: Optional[str] = None,
                 start: Optional[int] = None, end: Optional[int] = None,
                 strand : Optional[str] = None,
                 sequence: Optional[str] = None):
        self.id = exon_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.sequence = sequence or ""

        self.size = len(self.sequence) if self.sequence else 0

    def __repr__(self):
        return f"{self.id}:{self.chrom}\t{self.start}\t{self.end}\t{self.strand}"
    
    __str__ = __repr__

class CDS:
    """
    Represents an CDS in a gene.

    Attributes:
        id (str): The ID of the CDS.
        chrom (str): The chromosome where the CDS is located.
        start (int): The start position of the CDS.
        end (int): The end position of the CDS.
        strand (str): The strand of the CDS. "-" or "+"
        sequence (str): The sequence of the CDS, same with genome.
        size (int): The size of the CDS sequence.
        phase (int): The phase of the CDS.

    Methods:
        __repr__(): Returns a string representation of the CDS.
    """

    def __init__(self, cds_id: Optional[str] = None, 
                 chrom: Optional[str] = None,
                 start: Optional[int] = None, end: Optional[int] = None,
                 strand : Optional[str] = None,
                 sequence: Optional[str] = None,
                 phase: Optional[int] = None):
        self.id = cds_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.sequence = sequence or ""
        self._phase = phase

    def __repr__(self):
        return f"{self.id}:{self.chrom}\t{self.start}\t{self.end}\t{self.strand}\t{self.phase}"
    
    __str__ = __repr__

    @property
    def phase(self):
        return self._phase
    
    def set_phase(self, phase):
        self._phase = phase

class Orf:
    """
    Represents an Open Reading Frame (ORF) in gene annotation.

    Attributes:
        start (int): The start position of the ORF relative to the transcript.
        end (int): The end position of the ORF relative to the transcript.
        sequence (str): The sequence of the ORF, 5'->3'.
    """

    def __init__(self, start: Optional[int] = None, end: Optional[int] = None, sequence:Optional[str] = None, pep:Optional[str] = None):
        self.start = start
        self.end = end
        self.sequence = sequence or ""
        self._pep = pep
    
    @property
    def pep(self):
        if self._pep is None:
            self._pep = translate_cds(self.sequence)
        return self._pep
    
    def __str__(self):
        return f"{self.sequence}"
    
    __repr__ = __str__

    def __len__(self):
        return len(self.sequence)

    @classmethod
    def find_orf(cls, seq : str) :

        def find_start_codons(sequence):
            """Find all start codon positions in the given sequence."""
            return [i for i in range(0, len(sequence) - 2, 3) if sequence[i:i+3] == 'ATG']

        longest_orf: Optional[Orf] = None
        longest_len = 0

        seq = seq.upper()

        for frame in range(3):
            # Adjust length to be multiple of 3
            start = frame 
            end = len(seq) - (len(seq) - frame) % 3
            adjusted_seq = seq[start:end]
            trans = translate_cds(adjusted_seq)
            start_codons = find_start_codons(adjusted_seq)

            for start in start_codons:
                # Find the end of the ORF
                end = trans[start // 3:].find("*")
                if end != -1:
                    orf_len = (start + end * 3) - start
                    if orf_len > longest_len:
                        longest_len = orf_len
                        orf_start = start
                        orf_end = orf_start + orf_len + 3
   
                        longest_orf = cls(orf_start + frame, orf_end + frame , adjusted_seq[orf_start:orf_end] ,trans[orf_start // 3:orf_end // 3])

        return longest_orf

class Transcript:
    """
    Represents a transcript with its properties and methods.
    """

    def __init__(self, transcript_id: Optional[str] = None, 
                 chrom: Optional[str] = None,
                 start: Optional[int] = None,
                 end: Optional[int] = None,
                 strand: Optional[str] = None,
                 exon_list: Optional[List[Exon]] = None):
        """
        Initializes a Transcript object.

        Args:
            transcript_id (str, optional): The ID of the transcript. Defaults to None.
            chrom (str, optional): The chromosome of the transcript. Defaults to None.
            start (int, optional): The start position of the transcript. Defaults to None.
            end (int, optional): The end position of the transcript. Defaults to None.
            strand (str, optional): The strand of the transcript. Defaults to None.
            exon_list (List[Exon], optional): A list of exons associated with the transcript. Defaults to None.
        """
        self.id = transcript_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand

        self._sequence : str = None   
        self._exons = exon_list or []
        self._exon_num : int = len(self._exons)
        
        self._orf = None
        self._cds : List[CDS] = None
        

    def __repr__(self):
        """
        Returns a string representation of the Transcript object.

        Returns:
            str: The string representation of the Transcript object.
        """
        if self.chrom is not None:
            return f"{self.id}: {self.chrom}:{self.start}-{self.end}"
        else:
            return f"{self.id}"
         
    
    __str__ = __repr__

    @property
    def exons(self):
        """
        Returns the list of exons associated with the transcript.

        Returns:
            List[Exon]: The list of exons associated with the transcript.
        """
        return self._exons


    def add_exon(self, exon: Exon):
        """
        Adds an exon to the transcript.

        Args:
            exon (Exon): The exon to be added.
        """
        self._exons.append(exon)
        self._exon_num += 1

    @property
    def exon_num(self):
        """
        Returns the number of exons in the transcript.

        Returns:
            int: The number of exons in the transcript.
        """
        return self._exon_num
    
    def set_sequence(self):
        """
        Sets the sequence of the transcript by concatenating the sequences of its exons.
        """
        if not self._exons:
            self._sequence = ""
            return
        sorted_exons = sorted(self._exons, key=lambda exon: exon.start)
        sequence = "".join([exon.sequence for exon in sorted_exons ])

        if self.strand == "-":
            sequence = reverse_complement(sequence)

        self._sequence = sequence
    
    @property
    def sequence(self):
        """
        Returns the sequence of the transcript.

        Returns:
            str: The sequence of the transcript.
        """
        if self._sequence is None:
            self.set_sequence()
        return self._sequence
    
    @property
    def orf(self):
        """
        Returns the ORF of the transcript.

        Returns:
            Orf: The ORF of the transcript.
        """
        if self._orf is None:
            self._orf = Orf.find_orf(self.sequence)
        return self._orf
    
    
    @property
    def cds(self):
        """
        Returns the CDS of the transcript

        Returns:
            CDS: The CDS of the transcript.
        """
        if self._cds is not None:
            return self._cds
        
        cds_list = []

        orf_start = self.orf.start
        orf_end = self.orf.end
        # reverse the orf start and end if the transcript is on the negative strand
        if self.strand == "-":
            orf_start = len(self.sequence) - self.orf.end
            orf_end = len(self.sequence) - self.orf.start
    
        sorted_exons : List[Exon] = sorted(self._exons, key=lambda exon: exon.start)

        current_position = 0  # Position within the transcript
       
        for exon in sorted_exons:

            if current_position <= orf_start < current_position + exon.size:
                start_position = max(exon.start, exon.start + orf_start - current_position)
                end_position = min(exon.end, exon.start + orf_end - current_position)
                
                cds_sequence = exon.sequence[start_position-exon.start:end_position-exon.start]

                cds = CDS(f"cds_{self.id}", self.chrom, start_position, end_position, self.strand, cds_sequence)
                cds_list.append(cds)

                if orf_end <= current_position + exon.size:
                    break
            elif current_position < orf_end and orf_start < current_position:
                start_position = exon.start
                end_position = min(exon.end, exon.start + orf_end - current_position)
                
                cds_sequence = exon.sequence[start_position-exon.start:end_position-exon.start]
                
                cds = CDS(f"cds_{self.id}", self.chrom, start_position, end_position, self.strand, cds_sequence)
                cds_list.append(cds)

                if orf_end <= current_position + exon.size:
                    break
            current_position += exon.size

        # Filter out any CDS segments that may be outside of the ORF
        
        # set the phase of each cds
        cds_start = 0
        if self.strand == "-":
            cds_list = cds_list[::-1]
        for cds in cds_list:
            phase = ( 3- (cds_start % 3) ) % 3
            cds.set_phase(phase)
            cds_start += len(cds.sequence)
        
        self._cds = cds_list
        return self._cds
    
    def add_cds(self, cds: CDS):
        """
        Adds a CDS to the transcript.

        Args:
            cds (CDS): The CDS to be added.
        """
        if self._cds is None:
            self._cds = [ cds ]
        else:
            self._cds.append(cds)

class Gene:
    """
    Represents a gene with its attributes and associated transcripts.
    """

    def __init__(self, gene_id: Optional[str] = None, 
                 chrom: Optional[str] = None,
                 start: Optional[int] = None, end: Optional[int] = None,
                 strand: Optional[str] = None,
                 transcripts: Optional[List[Transcript] ]= None):
        """
        Initializes a Gene object.

        Args:
            gene_id (str, optional): The ID of the gene. Defaults to None.
            chrom (str, optional): The chromosome of the gene. Defaults to None.
            start (int, optional): The start position of the gene. Defaults to None.
            end (int, optional): The end position of the gene. Defaults to None.
            strand (str, optional): The strand of the gene. Defaults to None.
            transcripts (List[Transcript], optional): The list of associated transcripts. Defaults to None.
        """
        self.id = gene_id
        self._chrom = chrom
        self._start = start
        self._end = end
        self._strand = strand
        self._transcripts = transcripts
        self._sequence = None
        self._tx_num = 0

    def __repr__(self):
        
        if self._transcripts is None:
            return f"{self.id}: No transcripts"

        transcript_ids = [transcript.id for transcript in self._transcripts]
        all_tx = ",".join(transcript_ids)

        return f"{self.id}: { all_tx }"

    __str__ = __repr__

    @property
    def transcripts(self):
        """
        Returns the list of transcripts associated with the gene.

        Returns:
            List[Transcript]: The list of transcripts associated with the gene.
        """
        return self._transcripts

    def get_transcript(self, transcript_id: str) -> Optional[Transcript]:
        """
        Retrieves a transcript from the gene based on its ID.

        Args:
            transcript_id (str): The ID of the transcript to retrieve.

        Returns:
            Optional[Transcript]: The retrieved transcript, or None if not found.
        """
        if self._transcripts is None:
            return None

        for transcript in self._transcripts:
            if transcript.id == transcript_id:
                return transcript
        return None
    
    def add_transcript(self, transcript: Transcript) -> None:
        """
        Adds a transcript to the gene.
        
        Args:
            transcript (Transcript): The transcript to be added.
        """
        if self._transcripts is None:
            self._transcripts = [ transcript ]
            self._tx_num = 1
        else:
            self._transcripts.append(transcript)
            self._tx_num += 1

    
    def remove_transcript(self, transcript_id: str) -> bool:
        """
        Removes a transcript from the gene based on its ID.

        Args:
            transcript_id (str): The ID of the transcript to be removed.
        """
        for transcript in self._transcripts:
            if transcript.id == transcript_id:
                self._transcripts.remove(transcript)
                self._tx_num -= 1
                return True
        return False
    
    @property
    def chrom(self):
        if self._chrom is None:
            # get the chrom of the first transcript of the gene
            self._chrom = self._transcripts[0].chrom

        return self._chrom

    @property
    def start(self):
        if self._start is None:
            # get the start position of the first transcript of the gene
            all_start = [tx.start for tx in self._transcripts]
            # get the minimum start position
            self._start = min(all_start)

        return self._start
    
    @property
    def end(self):
        if self._end is None:
            # get the start position of the first transcript of the gene
            all_end = [tx.end for tx in self._transcripts]
            # get the minimum start position
            self._end = max(all_end)

        return self._end
    
    @property
    def strand(self):
        if self._strand is None:
            # get the strand of the first transcript of the gene
            self._strand = self._transcripts[0].strand

        return self._strand

    @property
    def tx_num(self):
        return self._tx_num

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

def gene_to_gff3(gene: Gene, source: str = 'wgap') -> str:
    gff3_lines = []

    # Gene line
    gff3_lines.append(
        f"{gene.chrom}\t{source}\tgene\t{gene.start + 1}\t{gene.end}\t.\t{gene.strand}\t.\tID={gene.id}"
    )

    # Transcript lines
    for transcript in gene.transcripts:
        gff3_lines.append(
            f"{gene.chrom}\t{source}\tmRNA\t{transcript.start+ 1}\t{transcript.end}\t.\t{transcript.strand}\t.\tID={transcript.id};Parent={gene.id}"
        )

        # Exon lines
        if len(transcript.exons) > 0:
            for exon in transcript.exons:
                gff3_lines.append(
                    f"{gene.chrom}\t{source}\texon\t{exon.start+ 1}\t{exon.end}\t.\t{exon.strand}\t.\tID={exon.id};Parent={transcript.id}"
                )

        # CDS lines
        if len(transcript.cds) > 0:
            for cds in transcript.cds:
                if cds.phase is None:
                    phase = "."
                else:
                    phase = str(cds.phase)
                gff3_lines.append(
                    f"{gene.chrom}\t{source}\tCDS\t{cds.start+ 1}\t{cds.end}\t.\t{cds.strand}\t{phase}\tID={cds.id};Parent={transcript.id}"
                )

    return "\n".join(gff3_lines)
