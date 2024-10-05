import sys
sys.path.append('modules_4main_script')
from modules_4main_script.module_for_filter_fastq import length_bounds, filter
def filter_fastq(fastq_file, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0):
    result_filtering = {}
    for seq_name, seq, quality in fastq_file.items():
        result_bounds = length_bounds(seq)
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

        if ((result_bounds>gc_bounds_l and result_bounds<gc_bounds_u) and
                (result_quality>quality_threshold) and
                (length_seq>length_bounds_l and length_seq<length_bounds_u)):
            result_filtering[seq_name] = (seq, quality)
    return result_filtering

