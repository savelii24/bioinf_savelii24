import argparse
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import logging
import os 

logging.basicConfig(
    filename='bioinf_ut.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)

logging.info("Script started")

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
        cysteine_positions = [
            f"Cysteine (C) position: {i + 1}"
            for i, aa in enumerate(self.sequence)
            if aa.upper() == "C"
        ]
        if cysteine_positions:
            return "\n".join(cysteine_positions)
        else:
            return "Cysteine (C) not found"

    def description(self):
        return "Amino Acid Sequence"


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0,
):
    """
    :param input_fastq: путь к входному файлу FASTQ
    :param output_fastq: путь к выходному файлу FASTQ
    :param gc_bounds: диапазон для фильтрации по проценту нуклеотидов G и C
    :param length_bounds: диапазон для фильтрации по длине последовательности
    :param quality_threshold: порог качества
    :return: фильтрация последовательностей по заданным параметрам
    """
    try:
        if not os.path.exists(input_fastq):
            raise FileNotFoundError(f"Input file '{input_fastq}' not found.")

        result_filtering = []

        with open(input_fastq, "r") as input_handle:
            for record in SeqIO.parse(input_handle, "fastq"):
                seq = record.seq
                quality = record.letter_annotations["phred_quality"]
                result_bounds = gc_fraction(seq) * 100
                result_quality = sum(quality) / len(quality) if quality else 0
                length_seq = len(seq)

                gc_bounds_l, gc_bounds_u = gc_bounds if isinstance(gc_bounds, tuple) else (0, gc_bounds)
                length_bounds_l, length_bounds_u = length_bounds if isinstance(length_bounds, tuple) else (0, length_bounds)

                if (
                    gc_bounds_l <= result_bounds <= gc_bounds_u and
                    length_bounds_l <= length_seq <= length_bounds_u and
                    result_quality >= quality_threshold
                ):
                    result_filtering.append(record)

        with open(output_fastq, "w") as output_handle:
            SeqIO.write(result_filtering, output_handle, "fastq")
        
        logging.info(f"Filtering complete. Output written to {output_fastq}.")
    except FileNotFoundError as e:
        logging.error(f"Error: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="FASTQ processor")
    parser.add_argument('input_fastq', help="Path to input FASTQ file")
    parser.add_argument('output_fastq', help="Path to output FASTQ file")
    parser.add_argument('--gc_bounds', type=str, default="0,100", help="GC bounds as 'min,max'")
    parser.add_argument('--length_bounds', type=str, default="0,100000", help="Length bounds as 'min,max'")
    parser.add_argument('--quality_threshold', type=int, default=0, help="Quality threshold")
    args = parser.parse_args()

    gc_bounds = tuple(map(int, args.gc_bounds.split(',')))
    length_bounds = tuple(map(int, args.length_bounds.split(',')))

    logging.info(f"Started filtering {args.input_fastq} with GC bounds {gc_bounds} and length bounds {length_bounds}.")
    
    filter_fastq(
        args.input_fastq,
        args.output_fastq,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=args.quality_threshold
    )

    logging.info("Script finished successfully.")
