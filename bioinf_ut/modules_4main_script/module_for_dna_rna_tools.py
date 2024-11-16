def transcribe(*seqs):
    """
    :param seqs: DNA-sequence(-s)
    :return: transcribed DNA-sequence
    """
    result = []
    seqs_seq = seqs[:-1]
    for seq in seqs_seq:
        transcribe_1_rule = seq.maketrans(
            "AaTtCcGg", "AaUuCcGg"
        )  # create rule for transcribing molecule
        result.append(seq.translate(transcribe_1_rule))
    if len(result) == 1:
        return result[0]
    else:
        return list(result)


def reverse(*seqs):
    """
    :param seqs: DNA-sequence(-s)
    :return: reversed DNA-sequence(-s)
    """
    result = []
    seqs_seq = seqs[:-1]
    for seq in seqs_seq:
        result.append(seq[::-1])  # return reverse sequence
    if len(result) == 1:
        return result[0]
    else:
        return list(result)


def complement(*seqs):
    """
    :param seqs: DNA-sequence(-s)
    :return: complementary DNA-sequence(-s)
    """
    result = []
    seqs_seq = seqs[:-1]
    for seq in seqs_seq:
        transcribe_1_rule = seq.maketrans(
            "AaTtCcGgUu", "TtAaGgCcAa"
        )  # create rule for complementary sequence
        result.append(seq.translate(transcribe_1_rule))
    if len(result) == 1:
        return result[0]
    else:
        return list(result)


def reverse_complement(*seqs):
    """
    :param seqs: DNA-sequence(-s)
    :return: reversed and complementary DNA-sequence(-s)
    """
    result = []
    seqs_seq = seqs[:-1]
    for seq in seqs_seq:
        reverse_complement_rule = seq.maketrans(
            "AaTtCcGgUu", "TtAaGgCcAa"
        )  # create rule for complementary sequence
        result.append(seq.translate(reverse_complement_rule)[::-1])
    if len(result) == 1:
        return result[0]
    else:
        return list(result)
