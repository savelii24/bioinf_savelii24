def length_bounds(seq):
    count_A = seq.count('A')
    count_T = seq.count('T')
    count_C = seq.count('C')
    count_G = seq.count('G')
    percent_gc = ((count_C+count_G)/(count_A+count_T+count_C+count_G)*100)
    return percent_gc

def quality_threshold(quality):
     seq_len = len(seq)
     score = 0
     for varience in quality:
        q_var = ord(varience)
        score = score + q_var
     final_score = score/seq_len
     return final_score
