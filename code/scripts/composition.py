import numpy as np
from scipy.stats import entropy
import math
import string
import crutil as util


# Amino acids
forbidden = ['B', 'J', 'O', 'U', 'X', 'Z']
aas = []
for a in string.ascii_uppercase:
    if a not in forbidden:
        aas.append(a)

def get_aa_freqs(seq):
    """
    Arguments:
    seq: a sequences (string representing a protein, drawn from the 20 naturally-occuring amino acids)
    Returns:
    seq_freqs_sorted: an interable of the frequency of each amino acid in the sequence.
        The indices match the categories in the same index in 'aas'. (20 natural amino acids in alphabetical order)
    """
    # Set sequence to uppercase
    seq = seq.upper()
    # Strip all characters except the amino acids
    seq = ''.join(c for c in seq if c in aas)
    
    # Parse sequence
    seq = list(seq)
    seq_aas, seq_counts = np.unique(seq, return_counts=True)
    zero_freq_aas = np.setxor1d(aas, seq_aas)
    # Get frequencies
    all_seq_aas = np.append(seq_aas, zero_freq_aas)
    seq_counts = np.append(seq_counts, np.zeros(zero_freq_aas.size))
    seq_freqs = seq_counts / len(seq)

    seq_freqs_sorted = seq_freqs[all_seq_aas.argsort()]
    return(seq_freqs_sorted)



def get_aa_counts(seq):
    """
    Arguments:
    seq: a sequences (string representing a protein, drawn from the 20 naturally-occuring amino acids)
    Returns:
    seq_counts_sorted: an interable of the counts of each amino acid in the sequence.
        The indices match the categories in the same index in 'aas'. (20 natural amino acids in alphabetical order)
    """
    # Set sequence to uppercase
    seq = seq.upper()
    # Strip all characters except the amino acids
    seq = ''.join(c for c in seq if c in aas)
    
    # Parse sequence
    seq = list(seq)
    seq_aas, seq_counts = np.unique(seq, return_counts=True)
    zero_freq_aas = np.setxor1d(aas, seq_aas)
    # Get frequencies
    all_seq_aas = np.append(seq_aas, zero_freq_aas)
    seq_counts = np.append(seq_counts, np.zeros(zero_freq_aas.size))

    seq_counts_sorted = seq_counts[all_seq_aas.argsort()]
    return(seq_counts_sorted)


def get_character_freqs(seq, counts = False, sort_output = True, gaps = False):
    """
    Arguments:
    seq: a sequences (string representing a protein, drawn from the 20 naturally-occuring amino acids)
    counts: bool, whether to return the counts in each category or the frequencies. False (returns freqs) by default.
    sort_output: bool, whether to return the frequencies or counts in alphabetically sorted order. True by default.
    If false, the categories will also be returned.
    Returns:
    aas: if unsorted, an iterable of the categories (amino acids) represented for each sequence
    freqs: an interable of the frequency of each category in the sequence.
        The indices match the categories in the same index in 'aas'.
        If counts is True, this will be raw counts rather than frequency.
        If sort_output is True only the frequencies will be returned, the assumed order is the
        amino acids in alphabetical order.
    """
    # Amino acids
    forbidden = ['B', 'J', 'O', 'U', 'X', 'Z']
    aas = []
    for a in string.ascii_uppercase:
        if a not in forbidden:
            aas.append(a)
    if gaps:
        aas.append('-')
    
    all_seq_aas = []
    all_seq_freqs = []
    
    # Parse sequence
    seq = list(seq)
    seq_aas, seq_counts = np.unique(seq, return_counts=True)
    zero_freq_aas = np.setxor1d(aas, seq_aas)

    seq_aas = np.append(seq_aas, zero_freq_aas)
    seq_freqs = np.append(seq_counts, np.zeros(zero_freq_aas.size))
    if not counts:
        seq_freqs = seq_freqs / len(seq)

    if sort_output:
        seq_freqs_sorted = seq_freqs[seq_aas.argsort()]
        all_seq_freqs.append(seq_freqs_sorted)

    else:
        all_seq_aas.append(seq_aas)
        all_seq_freqs.append(seq_freqs)
    
    if sort_output:
        return(all_seq_freqs)
    else:
        return(all_seq_aas, all_seq_freqs)


def count_aa(seq, states):
    '''
    Count the number of each amino acid in a sequence.
    '''
    aa_count = {}
    for state in states:
        aa_count[state] = 0
    for aa in seq:
        aa_count[aa] += 1
    return aa_count



def get_composition(seq, states, pseudocount=0.0, prior=None):
    '''
    Computes the composition of a sequence (seq) in terms of its constituent parts (states, iterable)

    Inputs:
        sequence: a string of characters

    Returns: An array representing the compositions
    '''
    n = len(states)
    states_dict = util.get_states_dict(states)
    dict_comp = count_aa(seq, states)
     
    for k in states_dict:
        if not prior is None:
            dict_comp[k] = (dict_comp[k] + pseudocount * prior[states_dict[k]] * n) /\
                            (len(seq) + pseudocount * n)
        else:
            dict_comp[k] = (dict_comp[k] + pseudocount) / (len(seq) + pseudocount * n)

    # dictionary to vector
    comp_vector = util.dict_to_vector(dict_comp, states)

    assert (np.sum(comp_vector) <= 1 + 1e-6) and (np.sum(comp_vector) >= 1 - 1e-6)

    return comp_vector


def calculate_seq_complexity(seq, typ='k1', norm=False):
    """
    Given a sequence (seq; str), return the sequence complexity in terms of 'Complexity' (type='k1', default)
    or entropy (type='k2')
    NOTE that this will only consider the part of the sequence that is one of the 20 natural amino acids --
    all other character will be stripped.
    Entropy normalization (only valid for type='k2') can be to the yeast proteome frequencies (norm='yeast'),
         equiprobable amino acids (norm='uniform'), or to a modified highly-charged entropy which assumes an
         charged fraction of 0.5 split equally between the four charged amino acid types, with the remaining
         0.5 split equally between all other amino acids.
    Assumes natural log base (could switch to base 20 later)
    """
    # Set sequence to uppercase
    seq = seq.upper()
    # Strip all characters except the amino acids
    seq = ''.join(c for c in seq if c in aas)

    L = len(seq)
    yeast_H_lang = 2.8947654337060476
    uniform_H_lang = 2.3261019050474445
    charged_H_lang = entropy([0.125]*4 + [0.03125]*16)
    
    
    # Calculate complexity vector
    sj = get_aa_counts(seq)
    sj[::-1].sort()
    
    if typ == "k1":
        num = math.factorial(L)
        denom = 1
        for n in sj:
            denom *= math.factorial(n)
        omega = num / denom
        complexity = (1 / L) * np.log(omega)
        
    if typ == "k2":
        complexity = 0
        for n in sj:
            if n != 0:
                complexity += (n / L) * (np.log(n / L))
        complexity = -1*complexity
        
        if norm == 'yeast':
            complexity = complexity / yeast_H_lang
        elif norm == 'uniform':
            complexity = complexity / uniform_H_lang
        elif norm == 'charged':
            complexity = complexity / charged_H_lang
        else:
            raise ValueError("Normalization ('norm') can be 'yeast', 'charged', or 'uniform'")

    
    return(complexity)


def calculate_hydropathy(s, scale='kd_norm'):
    """
    Given a sequence (a string or Seq object), calculate the mean hydropathy according to the normalzied Kyte-Doolittle scale
    (scale='kd_norm', default), the raw Kyte-Doolittle scale (scale='kd') or an externally supplied scale
    (must be a length-20 dictionary with amino acids coded as one-letter and upper case).
    The function will strip all non-amino acid characters.
    """
    kd = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "E": -3.5,
                      "Q": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
                      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
                      "Y": -1.3, "V": 4.2}
    
    if scale == 'kd_norm':
        vals = list(kd.values())
        scaled_vals = [v-np.min(vals) for v in vals]

        kd_scale_norm = {}
        for k in kd.keys():
            kd_scale_norm[k] = (kd[k] - np.min(vals)) / np.max(scaled_vals)
        
        scale = kd_scale_norm

    elif scale == 'kd':
        scale = kd
    
    seq = s.upper()
    seq = ''.join([i for i in seq if i in aas])
    if len(seq) == 0:
        raise ValueError("Invalid sequence supplied: {}".format(s))
    
    hydro = []
    for aa in seq:
        hydro.append(scale[aa])
    return(np.mean(hydro))