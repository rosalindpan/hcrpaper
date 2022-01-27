#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 17:30:20 2020

@author: Triandafillou

Utility functions for highly charged sequences
"""

import numpy as np
import string

forbidden = ['B', 'J', 'O', 'U', 'X', 'Z']
aas = []
for a in string.ascii_uppercase:
    if a not in forbidden:
        aas.append(a)


def charge_asym(seq):
    """
    Returns the charge asymmetry  
    """
    n = len(seq)
    f_minus = (seq.count('E') + seq.count('D')) / n
    f_plus = (seq.count('K') + seq.count('R')) / n

    try:
        sigma = (f_plus - f_minus)**2 / (f_plus + f_minus)
    except ZeroDivisionError:
        sigma = np.nan
    return(sigma)


def delta(seq, g):
    """
    """
    overall_sigma = charge_asym(seq)
    N = len(seq)

    delta = 0

    for i in np.arange(N-g+1):
        segment = seq[i:i+g]
        sigma_i = charge_asym(segment)
        if not np.isnan(sigma_i):
            delta += (sigma_i - overall_sigma)**2

    delta = delta / (N-g+1)
    return(delta)


def kappa(seq, g1=5, g2=6):
    nrep = int(len(seq)/2)
    #print('E'*nrep+'K'*nrep)
    max_d1 = delta('E'*nrep+'K'*nrep, g1)
    max_d2 = delta('E'*nrep+'K'*nrep, g2)

    k1 = delta(seq, g1) / max_d1
    k2 = delta(seq, g2) / max_d2

    return(np.mean([k1, k2]))


def mutate_move(mutation_prob, insertion_prob, deletion_prob):
    """
    allowed moves: duplication, deletion, mutation.
    """
    assert mutation_prob + insertion_prob + deletion_prob == 1.
    
    t2 = mutation_prob+insertion_prob
    
    r = np.random.rand()
    
    if r < mutation_prob:
        move="mut"
    elif (r > mutation_prob) & (r < t2):
        move="ins"
    else:
        move="del"
        
    return(move)



def evolve_sequence(start_seq, steps, mp=0.95, ip=0.025, dp=0.025):
    """
    """
    seq = list(start_seq)
    for i in np.arange(steps):
        site = int(len(seq)*np.random.random())
        move = mutate_move(mp, ip, dp)
        if move == "mut":
            if seq[site] == "E":
                seq[site] = "K"
            else:
                seq[site] = "E"
        elif move == "ins":
            insertion = "E" if np.random.random() > 0.5 else "K"
            seq.insert(site, insertion)
        else:
            seq.pop(site)
        #print(''.join(seq))

    return(''.join(seq))


def generate_seq_by_kappa(length, target_kappa, filename=None, save=False, verbose=False):
    """
    Returns a sequence with a kappa value close to the target kappa (within 3
    decimal points)
    """
    assert (target_kappa > 0) & (target_kappa < 1)

    starting_EK = ''.join(np.random.choice(['E', 'K'], length))

    while np.round(kappa(starting_EK), decimals=3) != np.round(target_kappa, decimals=3):
        starting_EK = ''.join(np.random.choice(['E', 'K'], length))

    if verbose:
        print('Final Sequence: {}'.format(starting_EK))
        print('Kappa = {}'.format(kappa(starting_EK)))

    #np.save('start_seq_kappa0.15_2.npy', starting_EK)
    if save:
        if filename is None:
            filename = "EKsequence_kappa{}.npy".format(target_kappa)

        np.save(filename, starting_EK)
        #np.save('start_seq_kappa0.15_2.npy', starting_EK)
    return(starting_EK)


def fraction_charged(sequence):
    """
    Given a string of single-letter amino acids (sequence), return the fraction of charged residues.
    """
    sequence = sequence.upper()
    n = len(sequence)
    #return((sequence.count('E') + sequence.count('K') + sequence.count('R') + sequence.count('D')) / n)
    try:
        fc = (sequence.count('E') + sequence.count('K') + sequence.count('R') + sequence.count('D')) / n
    except ZeroDivisionError:
        return(np.nan)

    return(fc)



def delta2(seq, g, N=None):
    """
    """
    overall_sigma = charge_asym(seq)
    assert not np.isnan(overall_sigma)
    
    if N is None:
        N = len(seq)

    delta = 0

    for i in np.arange(len(seq)-g+1):
        segment = seq[i:i+g]
        sigma_i = charge_asym(segment)
        if not np.isnan(sigma_i):
            delta += (sigma_i - overall_sigma)**2

    delta = delta / (N-g+1)
    return(delta)


def kappa_local(seq, g1=5, g2=6, v=False):
    nrep = int(len(seq)/2)
    
    dms = seq.count('E') + seq.count('D')
    dps = seq.count('K') + seq.count('R')
    dxs = int((len(seq) - (dms + dps))/3)
    
    #print('A'*dxs + 'E'*dms + 'A'*dxs + 'K'*dps + 'A'*dxs)
    max_d1 = delta2('A'*dxs + 'E'*dms + 'A'*dxs + 'K'*dps + 'A'*dxs, g1)
    max_d2 = delta2('A'*dxs + 'E'*dms + 'A'*dxs + 'K'*dps + 'A'*dxs, g2)
    
    if v:
        print("Max delta:")
        print(max_d1)
    
    k1 = delta2(seq, g1) / max_d1
    k2 = delta2(seq, g2) / max_d2

    return(np.mean([k1, k2]))

def sigma(f_plus, f_minus):
    """
    Returns the charge asymmetry  
    """
    try:
        s = (f_plus - f_minus)**2 / (f_plus + f_minus)
    except ZeroDivisionError:
        s = np.nan
    return(s)

def make_seq_by_fs(f_plus, f_minus, n,
                   shuffle=False, pad_value = 'A'):
    """
    Given a fraction of positive and negative charged residues, return a model sequence of length n
    with those properties of fraction charged.
    """
    assert f_plus + f_minus <= 1, "Fraction of charged residues must be sensical"
    n_m = int(f_minus*n)
    n_p = int(f_plus*n)
    seq = 'E'*n_m + 'K'*n_p + 'A'*(n - n_m - n_p)
    
    if shuffle:
        seq = list(seq)
        np.random.shuffle(seq)
        return(''.join(seq))
    else:
        return(seq)



def evolve_sequence_noindels(start_seq, steps, aas = ['E', 'K', 'A'], probs = [0.5, 0.5, 0]):
    """
    """
    seq = list(start_seq)
    for i in np.arange(steps):
        site = int(len(seq)*np.random.random())
        seq[site] = np.random.choice(aas, p=probs)

    return(''.join(seq))



def f_plus(sequence):
    """
    """
    sequence = sequence.upper()
    n = len(sequence)
    return((sequence.count('K') + sequence.count('R')) / n)

def f_minus(sequence):
    """
    """
    sequence = sequence.upper()
    n = len(sequence)
    return((sequence.count('E') + sequence.count('D')) / n)

def NPCR(sequence, window=None, av=False):
    """
    Given a sequence (str), return the net charge per residue, either of the whole protein (default), or by residue
    with a given window length (left-justified).
    Will return the signed value(s), set av to 'True' to return the absolute value.
    """
    if not window:
        if av:
            return(abs(f_plus(sequence) - f_minus(sequence)))
        else:
            return(f_plus(sequence) - f_minus(sequence))
    else:
        npcrs = []
        for r in np.arange(len(sequence) - window - 1):
            segment = sequence[r:r+window]
            npcrs.append(f_plus(segment) - f_minus(segment))
        
        if av:
            return([abs(x) for x in npcrs])
        else:
            return(npcrs)
        
    return(None)



