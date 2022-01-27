import scipy as sp

def count_aas(seq, aa):
    count = 0
    for i in seq:
        if i == aa:
            count += 1
    return count


def get_f_pos(seq):
    num_pos = count_aas(seq, 'K') + count_aas(seq, 'R') + count_aas(seq, 'H')
    f_pos = num_pos / len(seq)
    #print(f_pos)
    return f_pos


def get_f_neg(seq):
    num_neg = count_aas(seq, 'D') + count_aas(seq, 'E')
    f_neg = num_neg / len(seq)
    #print(f_neg)
    return f_neg


def get_sigma(seq):
    f_pos = get_f_pos(seq)
    f_neg = get_f_neg(seq)
    if (f_pos + f_neg) != 0:
        sigma = (f_pos - f_neg)**2 / (f_pos + f_neg)
    else: #?
        sigma = 0
    return sigma


def get_delta(seq, g):
    sigma = get_sigma(seq)
    sum_dev = 0
    i = 0
    while (i + g < (len(seq)+1)):
        blob = seq[i: i+g]
        sigma_i = get_sigma(blob)
        sum_dev += (sigma_i - sigma)**2
        i += 1
    delta = sum_dev / i
    return delta


def get_delta_max_1(seq, g):
    num_aa = int(len(seq) / 2) # truncates
    max_seq = 'K' * num_aa + 'E' * num_aa
    delta_max = get_delta(max_seq, g)
    return delta_max


def get_delta_max_2(seq, g):
    num_pos = count_aas(seq, 'K') + count_aas(seq, 'R')
    num_neg = count_aas(seq, 'D') + count_aas(seq, 'E')            
    num_uncharged = len(seq) - num_pos - num_neg
    
    max_seq = (num_uncharged // 2) * 'A'
    max_seq += num_pos * 'K'
    max_seq += num_neg * 'E'
    max_seq += (num_uncharged - num_uncharged // 2) * 'A'

    delta_max = get_delta(max_seq, g)
    
    return delta_max


def get_kappa(seq, g, alt_norm=False):
    delta = get_delta(seq, g)
    #print("delta =", delta, "when g =", g)
    if alt_norm:
        delta_max = get_delta_max_2(seq, g)
    else:
        delta_max = get_delta_max_1(seq, g)

    kappa = delta / delta_max
    return kappa


def get_avg_kappa(seq, g1=5, g2=6, alt_norm=False):
    #print("length of sequence:", len(seq))
    kappa1 = get_kappa(seq, g1, alt_norm)
    kappa2 = get_kappa(seq, g2, alt_norm)
    avg_kappa = (kappa1 + kappa2) / 2
    return avg_kappa


