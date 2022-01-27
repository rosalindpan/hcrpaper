import csv
from Bio import Seq, SeqIO
from datetime import datetime
import crutil as util
import string


def get_index_from_mapping(mapping, index):
    mapped = mapping[index]
    return mapped


def get_fractional_charge(seq):
    num_charged = 0
    charged_aas = ['K', 'R', 'D', 'E']
    for (i, aa) in enumerate(seq):
        if aa in charged_aas:
            num_charged += 1
    f_charge = num_charged / len(seq)
    return f_charge


def make_weights(window_size):
    weights = []
    len_half = window_size // 2
    if (window_size % 2 == 0):
        for i in range(len_half):
            weights.append(i/(len_half-1))
        full_weights = weights + weights[::-1] 
    if (window_size % 2 == 1):
        for i in range(len_half):
            weights.append(i/(len_half))
        full_weights = weights + [1.0] + weights[::-1] 
    return full_weights


def get_fractional_charge_weighted(seq, include_his=False):
    weights = make_weights(len(seq))
    num_charged = 0
    if include_his:
        charged_aas = ['K', 'R', 'D', 'E', 'H']
    else:
        charged_aas = ['K', 'R', 'D', 'E']
    for (i, aa) in enumerate(seq):
        if aa in charged_aas:
            num_charged += 1 * weights[i]
    f_charge = num_charged / sum(weights)
    return f_charge


def rolling_fractional_charges_weighted(seq, window_size, include_his=False):
    i = 0
    f_charges = []
    while (i + window_size <= len(seq)):
        blob = seq[i: i+window_size]
        f_charge = get_fractional_charge_weighted(blob, include_his=include_his)
        f_charges.append(f_charge)
        i += 1
    return f_charges


def is_above_threshold(f_charges, threshold):
    above_t = []
    for i, f in enumerate(f_charges):
        if (f >= threshold):
            above_t.append((i,f))
    return above_t

def identify_charged_regions(f_charges, threshold, tolerance=0):
    '''
    Minor issues:
    1. Does not identify regions of length 1
    2. list_of_regions returned does not include dips
    '''
    above_t = is_above_threshold(f_charges, threshold)
    list_of_regions = []
    new_entry = True
    total_diffs = 0
    counter = -1
    for i in range(1, len(above_t)):
        total_diffs += above_t[i][0] - above_t[i-1][0] - 1
        if total_diffs <= tolerance:
            if new_entry:
                list_of_regions.append([])
                counter += 1
                list_of_regions[counter].append(above_t[i-1])
                new_entry = False
            list_of_regions[counter].append(above_t[i])
        else:
            new_entry = True
            total_diffs = 0
    return list_of_regions


def select_min_length(regions, min_length):
    new_list = []
    for i in regions:
        if len(i) >= min_length:
            new_list.append(i)
    return new_list


def find_charged_regions(sequence, window_size, charge_threshold, min_length, tolerance):
    #print(mappings)
    gapless_seq = util.remove_gaps(sequence)
    f_charges = rolling_fractional_charges_weighted(gapless_seq, window_size)
    #print(f_charges)
    charged_regions = identify_charged_regions(f_charges, charge_threshold, tolerance=tolerance)
    long_charged_regions = select_min_length(charged_regions, min_length)
    return long_charged_regions


def range_to_seq(fullseq, region):
    seq_gapped = fullseq[region[0]:region[1]+1]
    seq = util.remove_gaps(seq_gapped)
    return seq


def filter_regions(sequence, charged_regions, charge_threshold):
    #print(sequence)
    result = []
    for i in charged_regions:
        seq = range_to_seq(sequence, i)
        frac_charge = get_fractional_charge(seq)
        #print(i)
        #print(frac_charge)
        if frac_charge >= charge_threshold:
            result.append(i)
    return result


def find_start_end(sequence, window_size, charge_threshold, min_length, tolerance,
					filter=True):
    sequence = str(sequence).replace("*", "")
    mappings = util.get_mapping(sequence)
    long_charged_regions = find_charged_regions(sequence, window_size,
                                charge_threshold, min_length, tolerance)
    rv = []
    for region in long_charged_regions:
        start = get_index_from_mapping(mappings, region[0][0])
        end = get_index_from_mapping(mappings, (region[-1][0] + window_size - 1))
        if (end == -1):
            end = get_index_from_mapping(mappings, len(sequence)-1)
        rv.append((start,end))

    #print(rv)
    # check if the fractional charge of regions detected are above a threshold
    if filter:
    	rv = filter_regions(sequence, rv, charge_threshold)

    # combine overlapping regions
    combined_regions = util.combine_regions(rv)

    # check again if the fractional charge of regions detected are above a threshold
    if filter:
    	combined_regions = filter_regions(sequence, combined_regions, charge_threshold)

    return combined_regions


def find_charged_seqs(sequence, window_size, charge_threshold, min_length, tolerance):
    charged_regions = find_start_end(sequence, window_size, charge_threshold,
                        min_length, tolerance)
    charged_seqs = []
    for i in charged_regions:
        charged_seq = range_to_seq(sequence, i)
        charged_seqs.append(charged_seq)
    return charged_seqs



# Helper functions for getting regions with high fraction of a null property from the proteome

## v0: letters from the first half of the alphabet

def find_alphabet_start_end(sequence, window_size, threshold, min_length, tolerance,
					filter=True):
    sequence = str(sequence).replace("*", "")
    mappings = util.get_mapping(sequence)
    long_regions = find_alphabet_regions(sequence, window_size,
                                threshold, min_length, tolerance)
    rv = []
    for region in long_regions:
        start = get_index_from_mapping(mappings, region[0][0])
        end = get_index_from_mapping(mappings, (region[-1][0] + window_size - 1))
        if (end == -1):
            end = get_index_from_mapping(mappings, len(sequence)-1)
        rv.append((start,end))

    #print(rv)
    # check if the fraction of letters in the first half of the alphabet of regions detected are above a threshold
    if filter:
    	rv = filter_alphabet_regions(sequence, rv, threshold)

    # combine overlapping regions
    combined_regions = util.combine_regions(rv)

    # check again if the fractional charge of regions detected are above a threshold
    if filter:
    	combined_regions = filter_alphabet_regions(sequence, combined_regions, threshold)

    return combined_regions


def find_alphabet_regions(sequence, window_size, threshold, min_length, tolerance):
    #print(mappings)
    gapless_seq = util.remove_gaps(sequence)
    f_firsts = rolling_alphabet_fraction_weighted(gapless_seq, window_size)
    #print(f_charges)
    regions = identify_charged_regions(f_firsts, threshold, tolerance=tolerance)
    long_regions = select_min_length(regions, min_length)
    return long_regions


def filter_alphabet_regions(sequence, regions, threshold):
    result = []
    for i in regions:
        seq = range_to_seq(sequence, i)
        frac_first = get_alphabet_fraction(seq)
        #print(i)
        #print(frac_charge)
        if frac_first >= threshold:
            result.append(i)
    return result


def rolling_alphabet_fraction_weighted(seq, window_size):
    i = 0
    f_firsts = []
    while (i + window_size <= len(seq)):
        blob = seq[i: i+window_size]
        f_first = get_alphabet_fraction_weighted(blob)
        f_firsts.append(f_first)
        i += 1
    return f_firsts


def get_alphabet_fraction_weighted(seq):
    weights = make_weights(len(seq))
    num_first_half = 0
    for (i, aa) in enumerate(seq.upper()):
        if aa in string.ascii_uppercase[0:13]:
            num_first_half += 1 * weights[i]
    f_first = num_first_half / sum(weights)
    return f_first

def get_alphabet_fraction(seq):
    num_first = 0
    for (i, aa) in enumerate(seq.upper()):
        if aa in string.ascii_uppercase[0:13]:
            num_first += 1
    f_first = num_first / len(seq)
    return f_first



## v1: enrichment for a specific subset of aas

def get_enriched_fraction(seq, target_aas):
    """
    Inputs:
    seq (str): a string of single-letter amino acids (can be upper or lower case).
    target_aas (list): the amino acids (upper or lower case) for which to get the fraction containing
    """
    num_targets = 0
    for (i, aa) in enumerate(seq.upper()):
        if aa in target_aas:
            num_targets += 1
    f_targets = num_targets / len(seq)
    return f_targets

def get_enriched_fraction_weighted(seq, target_aas):
    """
    Inputs:
    seq (str): a string of single-letter amino acids (can be upper or lower case).
    target_aas (list): the amino acids (upper or lower case) for which to get the fraction containing
    """
    weights = make_weights(len(seq))
    num_targets = 0
    for (i, aa) in enumerate(seq.upper()):
        if aa in target_aas:
            num_targets += 1 * weights[i]
    f_targets = num_targets / sum(weights)
    return f_targets


def rolling_enriched_fraction_weighted(seq, targets, window_size):
    i = 0
    f_targets = []
    while (i + window_size <= len(seq)):
        blob = seq[i: i+window_size]
        f_target = get_enriched_fraction_weighted(blob, targets)
        f_targets.append(f_target)
        i += 1
    return f_targets


def filter_enriched_regions(sequence, targets, regions, threshold):
    result = []
    for i in regions:
        seq = range_to_seq(sequence, i)
        frac_enriched = get_enriched_fraction(seq, targets)
        if frac_enriched >= threshold:
            result.append(i)
    return result


def find_enriched_regions(sequence, targets, window_size, threshold, min_length, tolerance):
    gapless_seq = util.remove_gaps(sequence)
    f_enriched = rolling_enriched_fraction_weighted(gapless_seq, targets, window_size)
    regions = identify_charged_regions(f_enriched, threshold, tolerance=tolerance)
    long_regions = select_min_length(regions, min_length)
    return long_regions



def find_enriched_start_end(sequence, targets, window_size, threshold, min_length, tolerance,
					filter=True):
    sequence = str(sequence).replace("*", "")
    mappings = util.get_mapping(sequence)
    long_regions = find_enriched_regions(sequence, targets, window_size,
                                threshold, min_length, tolerance)
    rv = []
    for region in long_regions:
        start = get_index_from_mapping(mappings, region[0][0])
        end = get_index_from_mapping(mappings, (region[-1][0] + window_size - 1))
        if (end == -1):
            end = get_index_from_mapping(mappings, len(sequence)-1)
        rv.append((start,end))

    #print(rv)
    # check if the fraction of enriched aas detected are above a threshold
    if filter:
    	rv = filter_enriched_regions(sequence, targets, rv, threshold)

    # combine overlapping regions
    combined_regions = util.combine_regions(rv)

    # check again if the fractional charge of regions detected are above a threshold
    if filter:
    	combined_regions = filter_enriched_regions(sequence, targets, combined_regions, threshold)

    return combined_regions