from modules_4main_script import module_for_dna_rna_tools

def run_dna_rna_tools(*seqs):
    operation = seqs[-1]
    seqs_seq = seqs[:-1]
    for seq in seqs_seq:
        if ((seq.count("t") + seq.count("T")) > 0) and (
                (seq.count("u") + seq.count("U")) > 0
        ):  # check the presence of uracil and thymine in sequence
            print(seq, "- урацил и тимин в одной молекуле")
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