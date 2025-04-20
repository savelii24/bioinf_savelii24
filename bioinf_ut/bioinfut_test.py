import unittest
from bioinf_ut import DNASequence, RNASequence, AminoAcidSequence, filter_fastq
from io import StringIO
import logging
import os

input_fastq = "example_fastq.fastq"
output_fastq = "test_output.fastq"

class TestBiologicalSequence(unittest.TestCase):
    #1
    def test_dna_sequence_creation(self):
        seq = "ATGCcccggaT"
        dna = DNASequence(seq)
        self.assertEqual(str(dna), seq)
        self.assertEqual(dna.description(), "DNA Sequence")
    #2
    def test_invalid_dna_sequence(self):
        seq = "ATGXqqwledj"
        with self.assertRaises(ValueError):
            DNASequence(seq)

class TestRNASequence(unittest.TestCase):
    #3
    def test_rna_sequence_creation(self):
        seq = "AUGC"
        rna = RNASequence(seq)
        self.assertEqual(str(rna), seq)
        self.assertEqual(rna.description(), "RNA Sequence")
    def test_invalid_rna_sequence(self):
        seq = "ATGX"
        with self.assertRaises(ValueError):
            RNASequence(seq)

class TestAminoAcidSequence(unittest.TestCase):
    #4
    def test_amino_acid_sequence_creation(self):
        seq = "ACDE"
        aa = AminoAcidSequence(seq)
        self.assertEqual(str(aa), seq)
        self.assertEqual(aa.description(), "Amino Acid Sequence")
    #5
    def test_find_cysteine_in_amino_acid_sequence(self):
        seq = "ACGAC"
        aa = AminoAcidSequence(seq)
        self.assertIn("Cysteine", aa.C_found())

class TestFilterFastq(unittest.TestCase):
    #6
    def test_filter_fastq_valid(self):
        gc_bounds = (30, 80)
        length_bounds = (100, 1000)
        quality_threshold = 30
        self.assertTrue(os.path.exists(input_fastq), f"File {input_fastq} not found")
        filter_fastq(input_fastq, output_fastq, gc_bounds, length_bounds, quality_threshold)
        self.assertTrue(os.path.exists(output_fastq))
        os.remove(output_fastq)
    
    #7
    def test_filter_fastq_invalid_file(self):
        gc_bounds = (30, 80)
        length_bounds = (100, 1000)
        quality_threshold = 30
        input_fastq = "non_existent_file.fastq"
        output_fastq = "output.fastq"
        with self.assertRaises(FileNotFoundError):
            filter_fastq(input_fastq, output_fastq, gc_bounds, length_bounds, quality_threshold)


class TestLogging(unittest.TestCase):
    def test_logging_info(self):
        log_stream = StringIO()
        logger = logging.getLogger("test_logger")
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler(log_stream)
        logger.addHandler(handler)
        logger.info("Test info log")
        handler.flush()
        log_output = log_stream.getvalue()
        self.assertIn("Test info log", log_output)
        logger.removeHandler(handler)

if __name__ == "__main__":
    unittest.main()
