input = '/bioinf_savelii24/bioinf_ut/example_files/example_fastq.fastq'
def read_input_fastq(input):
    with open(input) as file:
        input_fastq = file.readlines()
    return input_fastq
def fastq_filter(input_fastq):
    '''
    :param file: read file from your directory in FASTQ format
    :return: dictionary, where seq_name is a key, quality and seq is a value
    '''
    output_fastq_dict = {}
    for i in range(0, len(input_fastq), 4):
        seq_name = input_fastq[i].strip()
        seq = input_fastq[i + 1].strip()
        quality = input_fastq[i + 3].strip()
        output_fastq_dict[seq_name] = (seq, quality)
    return output_fastq_dict
output_path = '/bioinf_savelii24/bioinf_ut/filtered/result'
def save_to_file(output_path, output_fastq_dict):
    with open('/bioinf_savelii24/bioinf_ut/filtered/result', 'w') as file:
        for key, value in output_fastq_dict.items():
            file.write(f"{key}\n{value[0]}\n{value[1]}\n")
