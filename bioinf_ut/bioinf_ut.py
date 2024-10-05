from bioinf_ut.modules_4main_script.module_for_filter_fastq import (
    length_bounds_function,
)
from bioinf_ut.modules_4main_script.module_for_filter_fastq import filter


def filter_fastq(
    fastq_file, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0
):
    """
    :param fastq_file: fastq_file is dictionary with DNA-reads where key is sequence name (seq_name)
    and value is nucleotide sequence (seq) and quality of sequence in ASCII-code format (quality)
    :param gc_bounds: interval for filtering sequences by percent of C- and G- nucleotide in sequence
    :param length_bounds: value for filtering sequences by length (in nucleotide)
    :param quality_threshold: interval for filtering sequences by value of quality
    :return: filtering sequences by three parameters (gc_bounds, length_bounds and quality_threshold)
    """
    result_filtering = {}
    for seq_name, seq_1 in fastq_file.items():
        seq, quality = seq_1
        result_bounds = length_bounds_function(seq)
        result_quality = filter(quality)
        length_seq = len(seq)
        if isinstance(gc_bounds, int):
            gc_bounds_l, gc_bounds_u = 0, gc_bounds
        else:
            gc_bounds_l, gc_bounds_u = gc_bounds
        if isinstance(length_bounds, int):
            length_bounds_l, length_bounds_u = 0, length_bounds
        else:
            length_bounds_l, length_bounds_u = length_bounds

        if (
            gc_bounds_l <= result_bounds <= gc_bounds_u
            and gc_bounds_l <= length_seq <= length_bounds_u
            and result_quality >= quality_threshold
        ):
            result_filtering[seq_name] = (seq, quality)
    return result_filtering


from bioinf_ut.modules_4main_script.module_for_dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
)


def run_dna_rna_tools(*seqs):
    """
    :param seqs: DNA-sequence(-s)
    :return: DNA-sequence(-s) modified by one of four operation: transcribe, reverse, complement, reverse_complement
    function can check the presence of uracil and thymine in one sequence and return Error-message
    """
    operation = seqs[-1]
    seqs_seq = seqs[:-1]
    for seq in seqs_seq:
        if ((seq.count("t") + seq.count("T")) > 0) and (
            (seq.count("u") + seq.count("U")) > 0
        ):  # check the presence of uracil and thymine in one sequence
            print(seq, "- uracil and thymine in one sequence")
    for seq in seqs:
        if operation == "transcribe":
            final_result = transcribe(*seqs)
        elif operation == "reverse":
            final_result = reverse(*seqs)
        elif operation == "complement":
            final_result = complement(*seqs)
        elif operation == "reverse_complement":
            final_result = reverse_complement(*seqs)
    return final_result
