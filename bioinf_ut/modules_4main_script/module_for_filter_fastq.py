def length_bounds_function(seq):
    """
    :param seq: is DNA sequence from 'fastq_file'
    :return: percent of C- and G- nucleotide in seq
    """
    count_A = seq.count("A")
    count_T = seq.count("T")
    count_C = seq.count("C")
    count_G = seq.count("G")
    percent_gc = (count_C + count_G) / (count_A + count_T + count_C + count_G) * 100
    return percent_gc


def filter(quality):
    """
    :type quality: string
    :param quality: quality of DNA-read in
    ASCII-code format
    :return: value of quality in phred33 scale, border value is 0 (in a phred33 scale)
    """
    seq_len = len(quality)
    score = 0
    for var_1 in quality:
        q_var = ord(var_1) - 33
        score = score + q_var
    final_score = score / seq_len
    return final_score
