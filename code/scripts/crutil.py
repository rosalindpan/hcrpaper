'''
utility functions
'''

import numpy as np

def remove_gaps(seq):
    '''
    remove gaps ('-') in an amino acid sequence.
    input:
    seq - (str) an amino acid sequence
    output: 
    newseq - (str) the amino acid sequence with gaps removed
    '''
    return str(seq).replace('-','')


def get_mapping(seq):
    '''
    helper function for get_region_range; 
    maps indices of amino acid in the linear sequence to the indices of corresponding
    amino acids in a given alignment.
    input: 
    seq - (str) an amino acid sequence extracted directly from the alignment
    output: 
    mappings - (dict) one-to-one mapping between the corresponding indices
    '''
    mappings = {}
    counter = 0
    for (i,aa) in enumerate(seq):
        if aa != '-':
            mappings[counter] = i
            counter += 1
    return mappings


def combine_regions(regions):
    '''
    combines overlapping regions
    input:
    regions - (list of int tuples) list of regions
    output:
    combined_regiosn - (list of int tuples)
    '''
    regions = sorted(regions)
    # sorting is necessary: the rest of the function does not work if list is not sorted

    last_combined = regions
    combined_regions = []
    i = 0
    while i < len(regions):
        start = regions[i][0]
        tmp_end = regions[i][1]
        tmp_index = i
        for j in range(i+1, len(regions)):
            if regions[j][0] <= tmp_end:
                tmp_end = regions[j][1]
                tmp_index = j
        combined_regions.append((start, tmp_end))
        i = tmp_index + 1

    return combined_regions


def dict_to_vector(input_dict, states):
    '''
    Convert a dictionary to a 1d array in the order of the states.
    '''
    rv = []
    for state in states:
        if state in input_dict.keys():
            rv.append(input_dict[state])
        else:
            rv.append(0)
    vector = np.asarray(rv)
    return vector


def get_states_dict(states):
    '''
    Find a dictionary that maps each state in the state space to its index

    Analysis of real sequences uses all_aa as states
    '''
    states_dict = {}
    for i, state in enumerate(states):
        states_dict[state] = i

    return states_dict