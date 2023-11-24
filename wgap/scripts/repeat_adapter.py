from typing import Optional, Dict, List

class Repeat:
    """
    Represents a repetitive element identified by RepeatMasker in a DNA sequence.

    Attributes:
        id (str): The ID of the repeat.
        chrom (str): The chromosome where the repeat is located.
        start (int): The start position of the repeat.
        end (int): The end position of the repeat.
        strand (str): The strand of the repeat. "-" or "+"
        repeat_class (str): The class of the repeat, e.g., LINE, SINE, etc.
        repeat_family (str): The family of the repeat within its class.
        score (int): The score assigned by RepeatMasker, indicating the quality of the match.
        sequence (str): The sequence of the repeat, same with genome.
        size (int): The size of the repeat sequence.

    Methods:
        __repr__(): Returns a string representation of the repeat.
    """

    def __init__(self, repeat_id: Optional[str] = None, 
                 chrom: Optional[str] = None,
                 start: Optional[int] = None, end: Optional[int] = None,
                 strand: Optional[str] = None,
                 repeat_class: Optional[str] = None,
                 repeat_family: Optional[str] = None,
                 score: Optional[int] = None,
                 sequence: Optional[str] = None):
        self.id = repeat_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.repeat_class = repeat_class
        self.repeat_family = repeat_family
        self.score = score
        self.sequence = sequence

        self.size = len(self.sequence) if self.sequence else 0

    def __repr__(self):
        return (f"Repeat(ID: {self.id}, Chromosome: {self.chrom}, Start: {self.start}, "
                f"End: {self.end}, Strand: {self.strand}, Class: {self.repeat_class}, "
                f"Family: {self.repeat_family}, Score: {self.score}, Size: {self.size})")

    __str__ = __repr__
    

def repeat_masker_loader(data: List[str]) -> List[Repeat]:
    """
    Processes a list of strings representing the output from RepeatMasker.

    Args:
        data (List[str]): The RepeatMasker output lines.

    Returns:
        List[Repeat]: A list of Repeat objects parsed from the output.
    """
    
    repeats = []
    for line in data:
        # Skip header and empty lines
        if line.startswith('   SW') or line.strip() == '':
            continue

        # Extract fields from the line
        fields = line.split()

        # Parse fields and create Repeat object
        repeat_id = fields[15]
        chrom = fields[4]
        start = int(fields[5])
        end = int(fields[6])
        left = fields[7]
        strand = '+' if fields[8] == '+' else '-'
        repeat_sequence = fields[9]
        repeat_class_family = fields[10].split('/')
        repeat_class = repeat_class_family[0]
        repeat_family = repeat_class_family[1] if len(repeat_class_family) > 1 else None
        score = int(fields[0])
        sequence = None  # Sequence is not provided in the given data

        repeat = Repeat(repeat_id, chrom, start, end, strand, repeat_class, repeat_family, score, sequence)
        repeats.append(repeat)

    return repeats

# Sample data
data = [
    '   SW   perc perc perc  query         position in query           matching             repeat                position in repeat\n',
    'score   div. del. ins.  sequence      begin   end        (left)   repeat               class/family      begin   end    (left)     ID\n',
    '\n',
    '   14   13.4  2.9  2.9  ptg000006l_1    10305   10339 (1870736) + (AATT)n              Simple_repeat           1     35     (0)     1  \n',
    '   25   15.2  2.0  0.0  ptg000006l_1    10953   11003 (1870072) + (AGAATC)n            Simple_repeat           1     52     (0)     2  \n'
]

# Process the data
repeats = repeat_masker_loader(data)

# Display the results
for repeat in repeats:
    print(repeat)
