def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str):
    seqs = {}
    key = None
    seq = []
    with open(input_fasta, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line_1 = line.strip()
            if line_1.startswith('>'):
                if key:
                    seqs[key] = ''.join(seq)
                key = line
                seq = []
            else:
                seq.append(line_1)
        if key: seqs[key] = ''.join(seq)
    with open(output_fasta, "w") as outfile:
        for key, seq in seqs.items():
            outfile.write(f"{key}\n{seq}\n")
