from typing import Optional, Dict, List
from io import StringIO



SeqRecordDict = Dict[str, str]

# def read_fasta(fasta_file: str) -> SeqRecordDict:
#     fasta_dict = {}
#     for seq_record in SeqIO.parse(fasta_file, "fasta"):
#         fasta_dict[seq_record.id] = str(seq_record.seq)
#         return fasta_dict

def read_fasta(fasta_file: str) -> SeqRecordDict:
    fasta_dict = {}
    seq_builder = StringIO()  # 创建一个StringIO对象来有效地构建字符串
    seq_id = None
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                if seq_id is not None:
                    fasta_dict[seq_id] = seq_builder.getvalue()  # 保存之前的序列
                seq_id = line.strip()[1:]
                seq_builder = StringIO()  # 重置StringIO对象
            else:
                seq_builder.write(line.strip())  # 向StringIO对象中添加字符串
        if seq_id is not None:  # 不要忘记保存最后一个序列
            fasta_dict[seq_id] = seq_builder.getvalue()
    return fasta_dict

def write_fasta(fasta_dict: SeqRecordDict, output_file: str) -> None:
    with open(output_file, "w") as f:
        for seq_id, sequence in fasta_dict.items():
            f.write(f">{seq_id}\n")
            f.write(f"{sequence}\n")

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


