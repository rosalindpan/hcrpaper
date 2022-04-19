import numpy as np
import pandas as pd
import glob
import re
import gzip
import string
import os.path
from os import path
import mdtraj as md
import collections


"""
Utility functions for parsing and using alphafold outputs
"""

forbidden = ['B', 'J', 'O', 'U', 'X', 'Z']
aas = []
for a in string.ascii_uppercase:
    if a not in forbidden:
        aas.append(a)


def get_alphafold_fp(uniprot_id, path_to_af_data):
    """
    Given a uniprot reference number, search a directory (path_to_af_data) with the alphafold
    data, and return the full filepath to pdb file that predicts the structure of this orf.
    If multiple files are present only the first will be returned.
    """
    fpath = path_to_af_data + 'AF-' + str(uniprot_id) + '-F1-model_v2.pdb'

    return(fpath)


def is_gz_file(filepath):
	"""
	This function was sourced from Stack Overflow:
	https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
	"""
	with open(filepath, 'rb') as test_f:
		return test_f.read(2) == b'\x1f\x8b'


def read_bfactor_from_pdb(filepath):
    """
    Given a file path to a file in standard PDB format (filepath), return
    an array of the values in the b-factor column for each *residue* (NOT ATOM).
    This version of the function can handle both regular and gzipped files.
    Any non-ATOM lines are skipped.
    """
    # Test if file is gzipped
    if is_gz_file(filepath):
        with gzip.open(filepath) as file:
            resids = []
            pLDDTs = []
            for line in file:
                if line[0:6] == b"ATOM  ":
                    resnum = int(line[22:26].strip())
                    if resnum not in resids:
                        resids.append(resnum)
                        pLDDTs.append(float(line[61:67].strip()))
                    else:
                        continue
        if len(pLDDTs) == 0:
            print("Warning: file {} failed to generate a bfactor array".format(filepath))
        return(np.array(pLDDTs))
    
    else:
        with open(filepath) as file:
            resids = []
            pLDDTs = []
            for line in file:
                if line[0:6] == "ATOM  ":
                    resnum = int(line[22:26].strip())
                    if resnum not in resids:
                        resids.append(resnum)
                        pLDDTs.append(float(line[61:67].strip()))
                    else:
                        continue
        if len(pLDDTs) == 0:
            print("Warning: file {} failed to generate a bfactor array".format(filepath))        
        return(np.array(pLDDTs))


def read_af_output(fdir, uniprot_id):
    '''
    read in AlphaFold pdb file as a mdtraj file.
    input:
    fdir - (str) filepath to AlphaFold data
    uniprot_id - (str) 
    output:
    mdtraj pdb file
    '''
    fpath = fdir + 'AF-' + str(uniprot_id) + '-F1-model_v1.pdb'
    if path.exists(fpath):
        pdb = md.load(fpath)
    else:
        pdb = None
    return pdb



def read_seq_from_pdb(filepath):
    """
    Given a file path to a file in standard PDB format (filepath), return
    its amino acid sequence as a string of one letter amino acids
    This version of the function can handle both regular and gzipped files.
    """
    
    three_to_one_map =  {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', \
                         'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', \
                         'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', \
                         'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y', '-':'-'}
    # Test if file is gzipped
    if is_gz_file(filepath):
        with gzip.open(filepath) as file:
            resids = []
            pLDDTs = []
            for line in file:
                if line[0:6] == b"SEQRES":
                    resids+=line[19:71].decode("utf-8").strip().split()
        
        return(''.join([three_to_one_map[aa] for aa in resids]))
    
    else:
        with open(filepath) as file:
            resids = []
            pLDDTs = []
            for line in file:
                if line[0:6] == "SEQRES":
                    resids+=line[19:71].strip().split()
        
        return(''.join([three_to_one_map[aa] for aa in resids]))


def read_af_output(fdir, uniprot_id):
    '''
    read in AlphaFold pdb file as a mdtraj file.
    input:
    fdir - (str) filepath to AlphaFold data
    uniprot_id - (str) 
    output:
    mdtraj pdb file
    '''
    fpath = fdir + 'AF-' + str(uniprot_id) + '-F1-model_v1.pdb'
    if path.exists(fpath):
        pdb = md.load(fpath)
    else:
        pdb = None
    return pdb


def get_percent_helix(fdir, uniprot_id, left_bound, right_bound):
    af_pdb = read_af_output(fdir, uniprot_id)
    if not af_pdb is None:
        ss = md.compute_dssp(af_pdb, simplified=True)[0]
        region_ss = ss[left_bound:(right_bound+1)]
        helix_cnt = collections.Counter(region_ss)['H']
        helix_p = helix_cnt / (right_bound - left_bound + 1)
    else:
        helix_p = None
    return helix_p


def get_freq_values_in_range(arr, lower_bound, upper_bound):
    '''
    helper function for get_disorder_label
    '''
    cnt_in_range = 0
    for i in arr:
        if (i >= lower_bound) and (i <= upper_bound):
            cnt_in_range += 1
    return cnt_in_range / len(arr) 


def get_structure_label(fdir, uniprot_id, left_bound, right_bound,
                        helical_cutoff=0.8, bfactor_cutoff=0.5):
    fpath = fdir + 'AF-' + str(uniprot_id) + '-F1-model_v1.pdb'
    bfactor = read_bfactor_from_pdb(fpath)[left_bound:(right_bound+1)]
    len_region = right_bound - left_bound + 1
    disorder_freq = get_freq_values_in_range(bfactor, 0, 50)
    order_freq = get_freq_values_in_range(bfactor, 70, 100)
    p_helix = get_percent_helix(fdir, uniprot_id, left_bound, right_bound)
    if disorder_freq >= bfactor_cutoff:
        label = 'disordered'
    elif (order_freq >= bfactor_cutoff) and (p_helix >= helical_cutoff):
        label = 'helix'
    else:
        label = 'unclassified'
    return label



### DELTE THIS AFTER TESTING ALL FIGURE NOTEBOOKS


# def get_percent_helix(ss, bfactor, len_region):
#     cnt_helix = 0
#     for i, label in enumerate(ss):
#         if (label == 'H') and (bfactor[i] >= 70):
#             cnt_helix += 1
#     return cnt_helix / len_region


# def get_percent_disorder(ss, bfactor, len_region):
#     cnt_disorder = 0
#     for i, label in enumerate(ss):
#         if ((label == 'C') and (bfactor[i] >= 70)) or (bfactor[i] < 70):
#             cnt_disorder += 1
#     return cnt_disorder / len_region


# def get_freq_values_in_range(arr, lower_bound, upper_bound):
#     '''
#     helper function for get_disorder_label
#     '''
#     cnt_in_range = 0
#     for i in arr:
#         if (i >= lower_bound) and (i <= upper_bound):
#             cnt_in_range += 1
#     return cnt_in_range / len(arr) 


# def get_structure_label(fdir, uniprot_id, left_bound, right_bound,
#                         helical_cutoff=0.6, disorder_cutoff=0.6):
#     fpath = fdir + 'AF-' + str(uniprot_id) + '-F1-model_v1.pdb'
#     if path.exists(fpath):
#         af_pdb = read_af_output(fdir, uniprot_id)
#         ss = md.compute_dssp(af_pdb, simplified=True)[0]
#         region_ss = ss[left_bound:(right_bound+1)]
#         bfactor = read_bfactor_from_pdb(fpath)[left_bound:(right_bound+1)]
#         len_region = right_bound - left_bound + 1

#         p_helix = get_percent_helix(region_ss, bfactor, len_region)
#         p_disorder = get_percent_disorder(region_ss, bfactor, len_region)

#         if p_disorder >= disorder_cutoff:
#             label = 'disordered'
#         elif p_helix >= helical_cutoff:
#             label = 'helix'
#         else:
#             label = 'unclassified'
#     else:
#         label = None
#     return label

# # below 70% or above 70% and predicted coil


