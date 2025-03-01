from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

class BiologicalSequence(ABC):
    def __init__(self, sequence: str, alphabet: set):
        self.sequence = sequence
        self.alphabet = alphabet
        if not self._is_valid():
            raise ValueError(f"Invalid sequence: {sequence}")
    def __len__(self):
        return len(self.sequence)
    def __getitem__(self, index):
        return self.sequence[index]
    def __str__(self):
        return self.sequence
    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"
    def _is_valid(self):
        return set(self.sequence).issubset(self.alphabet)
    @abstractmethod
    def description(self):
        pass

class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str, alphabet: set):
        super().__init__(sequence, alphabet)
    def complement(self):
        transcribe_rule = str.maketrans("AaTtCcGgUu", "TtAaGgCcAa") 
        return self.sequence.translate(transcribe_rule)
    def reverse(self):
        return self.sequence[::-1]
    def reverse_complement(self):
        transcribe_rule = str.maketrans("AaTtCcGgUu", "AaTtCcGgAa")
        return self.sequence.translate(transcribe_rule)[::-1]
    def description(self):
        return "Nucleic Acid Sequence"
    
class DNASequence(NucleicAcidSequence):
    DNA_ALPHABET = {"A", "T", "G", "C", "a", "t", "g", "c"}
    def __init__(self, sequence: str):
        super().__init__(sequence, self.DNA_ALPHABET)
    def transcribe(self):
        transcribe_rule = str.maketrans("AaTtCcGg", "AaUuCcGg")
        return self.sequence.translate(transcribe_rule)
    def description(self):
        return "DNA Sequence"

class RNASequence(NucleicAcidSequence):
    RNA_ALPHABET = {"A", "U", "G", "C", "a", "u", "g", "c"}
    def __init__(self, sequence: str):
        super().__init__(sequence, self.RNA_ALPHABET)
    def description(self):
        return "RNA Sequence"

class AminoAcidSequence(BiologicalSequence):
    AMINO_ACID_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")
    def __init__(self, sequence: str):
        super().__init__(sequence, self.AMINO_ACID_ALPHABET)
    def C_found(self):
        cysteine_positions = [f"Сysteine (С) posiotion: {i + 1}" for i, aa in enumerate(self.sequence) if aa.upper() == "C"]
        if cysteine_positions:
            return "\n".join(cysteine_positions)
        else:
            return "No cysteine (C) in transcript"
    def description(self):
        return "Amino Acid Sequence"


def filter_fastq(
    input_fastq: str,  output_fastq: str, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
    """
    :param fastq_file: fastq_file is dictionary with DNA-reads where key is sequence name (seq_name)
    and value is nucleotide sequence (seq) and quality of sequence in ASCII-code format (quality)
    :param gc_bounds: interval for filtering sequences by percent of C- and G- nucleotide in sequence
    :param length_bounds: value for filtering sequences by length (in nucleotide)
    :param quality_threshold: interval for filtering sequences by value of quality
    :return: filtering sequences by three parameters (gc_bounds, length_bounds and quality_threshold)
    """
    result_filtering = []
    for record in SeqIO.parse(input_fastq, "fastq"):
        seq = record.seq
        quality = record.letter_annotations["phred_quality"]
        result_bounds = gc_fraction(seq) * 100 
        result_quality = sum(quality) / len(quality) if len(quality) > 0 else 0
        length_seq = len(seq)
        if isinstance(gc_bounds, tuple):
            gc_bounds_l, gc_bounds_u = gc_bounds
        else:
            gc_bounds_l, gc_bounds_u = 0, gc_bounds
        if isinstance(length_bounds, tuple):
            length_bounds_l, length_bounds_u = length_bounds
        else:
            length_bounds_l, length_bounds_u = 0, length_bounds
        if (
            gc_bounds_l <= result_bounds <= gc_bounds_u
            and length_bounds_l <= length_seq <= length_bounds_u
            and result_quality >= quality_threshold
        ):
            result_filtering.append(record)
    SeqIO.write(result_filtering, output_fastq, "fastq")
