from typing import Optional, Dict, List
from Bio import SeqIO

SeqRecordDict = Dict[str, str]

def read_fasta(fasta_file: str) -> SeqRecordDict:
    fasta_dict = {}
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[seq_record.id] = str(seq_record.seq)
    return fasta_dict

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

## helpul functions
def reverse_complement(sequence: str) -> str:
    """
    Returns the reverse complement of a sequence.

    Args:
        sequence (str): The sequence to be reverse complemented.

    Returns:
        str: The reverse complement of the sequence.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C",
                  "a": "t", "t": "a", "c": "g", "g": "c"}
    
    try:
        return "".join([complement[base] for base in sequence[::-1]])
    except KeyError:
        print(f"Error: Invalid base in sequence: {sequence}")

# Default codon table (Standard Genetic Code)
STANDARD_CODON_TABLE = {
    "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
    "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
    "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
    "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
    "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
    "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
    "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
    "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
    "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
    "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
    "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
    "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
    "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
    "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
    "TAC":"Y", "TAT":"Y", "TAA":"*", "TAG":"*",
    "TGC":"C", "TGT":"C", "TGA":"*", "TGG":"W",
}

def translate_cds(dna_sequence, codon_table: Dict[str, str] = STANDARD_CODON_TABLE ) -> str:
    """
    Translate a DNA sequence into a protein sequence using a given codon table.

    Args:
    dna_sequence (str): The DNA sequence to be translated.
    codon_table (dict, optional): A dictionary representing the codon table. 
                                  If None, a default codon table is used.

    Returns:
    str: The translated protein sequence.
    """
    # Ensure the DNA sequence length is a multiple of 3
    if len(dna_sequence) % 3 != 0:
        raise ValueError("DNA sequence length is not a multiple of 3.")
    
    protein_sequence = ""
    # Translate each codon into an amino acid
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        protein_sequence += codon_table.get(codon, "?")  # '?' for unknown codons

    return protein_sequence

# Exon -> Transcript -> Gene

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
        self.sequence = sequence

        self.size = len(self.sequence)

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
        self.sequence = sequence
        self._phase = phase

    def __repr__(self):
        return f"{self.id}:{self.chrom}\t{self.start}\t{self.end}\t{self.strand}\t{self.phase}"
    
    __str__ = __repr__

    def set_phase(self, phase):
        self._phase = phase

    @property
    def phase(self):
        return self._phase

class Orf:
    """
    Represents an Open Reading Frame (ORF) in gene annotation.

    Attributes:
        start (int): The start position of the ORF relative to the transcript.
        end (int): The end position of the ORF relative to the transcript.
        sequence (str): The sequence of the ORF, 5'->3'.
    """

    def __init__(self, start: Optional[int] = None, end: Optional[int] = None, cds:Optional[str] = None, pep:Optional[str] = None):
        self.start = start
        self.end = end
        self.cds = cds
        self._pep = pep
    
    @property
    def pep(self):
        if self._pep is None:
            self._pep = translate_cds(self.cds)
        return self._pep
    
    def __str__(self):
        return f"{self.cds}"
    
    __repr__ = __str__

    def __len__(self):
        return len(self.cds)

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
                 start: Optional[int] = None, end: Optional[int] = None,
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
        self._exons = exon_list
        self._sequence : str = None
        self._exon_num : int = 0
        
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
        if self._exons is None:
            self._exons = [ exon ]
            self._exon_num = 1
        else:
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
        Returns the CDS of the transcript.

        Returns:
            CDS: The CDS of the transcript.
        """
        if self._cds is not None:
            return self._cds
        
        cds_list = []

        orf_start = self.orf.start
        orf_end = self.orf.end
        if self.strand == "-":
            orf_start = len(self.sequence) - self.orf.end
            orf_end = len(self.sequence) - self.orf.start

        sorted_exon : List[Exon] = sorted(self._exons, key=lambda exon: exon.start)
        exon_relative_start = 0
        exon_relative_end = 0
        for exon in sorted_exon:
            exon_relative_start = exon_relative_end
            exon_relative_end = exon_relative_start +  exon.size

            # compare each exon relative position with orf start and end
            if orf_start >= exon_relative_start and orf_start <= exon_relative_end:
                # orf start is in this exon
                cds_start = exon.start + orf_start
                cds_end = exon.end
                cds_sequence = exon.sequence[cds_start-exon.start:]
                cds = CDS(f"cds_{self.id}", self.chrom, cds_start, cds_end, self.strand, cds_sequence)
                cds_list.append(cds)
            elif orf_end >= exon_relative_start and orf_end <= exon_relative_end:
                # orf end is in this exon
                cds_start = exon.start
                cds_end = exon.start + (  orf_end - exon_relative_start )
                cds_sequence = exon.sequence[:cds_end-exon.start]
                cds = CDS(f"cds_{self.id}", self.chrom, cds_start, cds_end, self.strand, cds_sequence)
                cds_list.append(cds)
            elif orf_start <= exon_relative_start and orf_end >= exon_relative_end:
                # orf is in this exon
                cds_start = exon.start
                cds_end = exon.end
                cds_sequence = exon.sequence
                cds = CDS(f"cds_{self.id}", self.chrom, cds_start, cds_end, self.strand, cds_sequence)
                cds_list.append(cds)
        
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

# Loaders
def stringtie_loader(stringtie_gtf: str, fasta_dict: Optional[SeqRecordDict] = None) -> Dict[str, Gene]:
    
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

def gff3_loader(gff3_file: str, fasta_dict: Optional[SeqRecordDict] = None) -> Dict[str, Gene]:
    
    gene_dict : Dict[str, Gene] = {} 
    exon_dict : Dict[str, Exon] = {}
    
    for line in open(gff3_file):
        if line.startswith("#"):
            continue
        line = line.strip()

        #parts = re.split(r'\s+', line.strip())
        parts = line.split("\t")
        if len(parts) != 9:
            continue  # Ignore malformed lines

        chrom, _, feature_type, start, end, _, strand, _, attributes = parts
        # gff is 1-based, convert to 0-based
        start, end = int(start) - 1, int(end)

        attributes_dict = dict(item.split("=") for item in attributes.split(";") if "=" in item)

        if chrom not in fasta_dict:
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
        for exon in transcript.exons:
            gff3_lines.append(
                f"{gene.chrom}\t{source}\texon\t{exon.start+ 1}\t{exon.end}\t.\t{exon.strand}\t.\tID={exon.id};Parent={transcript.id}"
            )

        # CDS lines
        for cds in transcript.cds:
            gff3_lines.append(
                f"{gene.chrom}\t{source}\tCDS\t{cds.start+ 1}\t{cds.end}\t.\t{cds.strand}\t{cds.phase}\tID={cds.id};Parent={transcript.id}"
            )

    return "\n".join(gff3_lines)


def main(args):
    fasta = args.ref
    gene_dict = None

    if args.format == "stringtie":
        gene_dict = stringtie_loader(args.gff, fasta)

    min_exon_num = args.min_exon_num
    max_exon_num = args.max_exon_num
    min_orf_size = args.min_orf_size

    # step1: parse gene_dict to filter transcripts with exon number < 2 or exon number > 1000
    gene_del_list = []
    for gene in gene_dict.values():
        if gene.transcripts is None:
            gene_del_list.append(gene.id)
            continue
        for tx in gene.transcripts:
            if tx.exon_num < min_exon_num or tx.exon_num > max_exon_num:
                gene.remove_transcript(tx.id)
        if gene.tx_num == 0:
            gene_del_list.append(gene.id)

    # delete genes with no transcripts
    for gene_id in gene_del_list:
        del gene_dict[gene_id]

    # optional: filter genes which overlap with repeat regions
    if args.rm_out is not None:
        pass

    # step2: select the transcript with longest ORF for each gene and size > 300
    new_gene_dict :Dict[str, Gene] = {}
    for gene in gene_dict.values():
        longest_orf = None
        longest_orf_len = 0
        longest_tx = None
        for tx in gene.transcripts:
            if tx.orf is not None and len(tx.orf.cds) > longest_orf_len:
                longest_orf = tx.orf
                longest_orf_len = len(tx.orf.cds)
                longest_tx = tx
            
        if longest_orf is not None and len(longest_orf.cds) > min_orf_size :
            new_gene_dict[gene.id] = Gene(gene.id, gene.chrom, gene.start, gene.end, gene.strand)
            new_gene_dict[gene.id].add_transcript(longest_tx)

    # step3: export the new gene_dict to gff3 and fasta
    gff3_out = f"{args.prefix}.gff3"
    cds_out = f"{args.prefix}.cds.fa"
    with open(gff3_out, "w") as gff3_file:
        for gene in new_gene_dict.values():
            if gene.transcripts is None:
                continue
            gff3_file.write(gene_to_gff3(gene) + "\n")

    # step3: export the cds sequence to fasta
    with open(cds_out, "w") as cds_file:
        for gene in new_gene_dict.values():
            if gene.transcripts is None:
                continue
            for tx in gene.transcripts:
                cds_file.write(f">{tx.id}\n{tx.orf.cds}\n")

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description='gene annotation adapter')
    parser.add_argument('--rm_out', type=str, help='repeatmasker output file, for example, genome.fa.out')
    parser.add_argument('-m', '--min_exon_num', type=int, default=2, help='min exon number')
    parser.add_argument('-M', '--max_exon_num', type=int, default=1000, help='max exon number')
    parser.add_argument('-s', '--min_orf_size', type=int, default=300, help='min orf size')
    parser.add_argument('-p','--prefix', type=str, help='prefix of output file', default='wgap')
    parser.add_argument('-f', '--format', type=str, help='format of input gff file', default='stringtie')
    parser.add_argument('ref', type=str, help='fasta file')
    parser.add_argument('gff', type=str, help='input gtf file')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argparse()
    main(args)

