from typing import Optional, Dict, List
from .seq_utils import reverse_complement, translate_cds


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
    
    def to_gff3(self, source='wgap') -> str:
        # Validate required gene attributes
        required_attributes = ['chrom', 'start', 'end', 'strand', 'id', 'transcripts']
        for attr in required_attributes:
            if not hasattr(self, attr):
                raise ValueError(f"Gene object is missing required attribute: {attr}")

        # Initialize GFF3 lines list
        gff3_lines = []

        # Gene line
        gff3_lines.append(
            f"{self.chrom}\t{source}\tgene\t{self.start + 1}\t{self.end}\t.\t{self.strand}\t.\tID=gene:{self.id}"
        )

        # Validate and append transcript lines
        for transcript in self.transcripts:
            # Validate transcript attributes
            for attr in ['start', 'end', 'strand', 'id', 'exons' ]:
                if not hasattr(transcript, attr):
                    raise ValueError(f"Transcript object is missing required attribute: {attr}")

            gff3_lines.append(
                f"{self.chrom}\t{source}\tmRNA\t{transcript.start + 1}\t{transcript.end}\t.\t{transcript.strand}\t.\tID=transcript:{transcript.id};Parent=gene:{self.id}"
            )

            # Validate, format, and append exon lines
            if transcript.exons:
                for exon in transcript.exons:
                    gff3_lines.append(
                        f"{self.chrom}\t{source}\texon\t{exon.start + 1}\t{exon.end}\t.\t{exon.strand}\t.\tID=exon:{exon.id};Parent=transcript:{transcript.id}"
                    )

            # Validate, format, and append CDS lines
            if transcript.cds:
                for cds in transcript.cds:
                    phase = '.' if cds.phase is None else str(cds.phase)
                    gff3_lines.append(
                        f"{self.chrom}\t{source}\tCDS\t{cds.start + 1}\t{cds.end}\t.\t{cds.strand}\t{phase}\tID=cds:{cds.id};Parent=transcript:{transcript.id}"
                    )

        return "\n".join(gff3_lines)
