def read_input_fastq(input_path: str):
    """
    :param input_path: path to fili with fastq-reads
    :return: reading file
    """
    with open(input_path) as file:
        input_fastq = file.readlines()
    return input_fastq


def fastq_dict(input_fastq):
    """
    :param file: read file from your directory in FASTQ format
    :return: dictionary, where seq_name is a key, quality and seq is a value
    """
    output_fastq_dict = {}
    for i in range(0, len(input_fastq), 4):
        seq_name = input_fastq[i].strip()
        seq = input_fastq[i + 1].strip()
        quality = input_fastq[i + 3].strip()
        output_fastq_dict[seq_name] = (seq, quality)
    return output_fastq_dict


def save_fastq(output_path: str, output_fastq_dict: dict[str, tuple[str, str]]):
    """
    :param output_path: path to file with output
    :param output_fastq_dict: fastq reads in dictionary format
    :return: file with result
    """
    with open(output_path, "w") as file:
        for key, value in output_fastq_dict.items():
            file.write(f"{key}\n{value[0]}\n{value[1]}\n")
