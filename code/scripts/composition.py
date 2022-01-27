import numpy as np
from scipy.stats import entropy
import math
import string


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