import random
import numpy as np
from scipy import stats
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import crutil

'''
utility functions
'''

def load_msa(fdir, orf):
    '''
    read in msa as an alignment object.
    **NOTE: this function assumes specific formatting**
    input:
    fdir - (str) directory of msa file
    orf - (str) orf name
    output:
    msa - (AlignIO object)
    '''
    msa = AlignIO.read(open(fdir+str(orf)+'_aybrah.fa'), "fasta")
    return msa

'''
auxiliary functions
'''

def get_region_range(msa, left_bound, right_bound,
                     refseq=None):
    '''
    map a region range in reference sequence to range in msa.
    **NOTE: the function might not work well for non-aybrah alignment**
    input:
    msa - (alignIO object) msa
    left_bound - (int) left boundary of region in reference sequence
    right_bound - (int) right boundary of region in reference sequence
    refseq - (str, optional) reference sequence; does not have to be full seq
    output: 
    msa_left_bound - (int) left region boundary in msa
    msa_right_bound - (int) right region boundary in msa
    '''

    if refseq is not None: 
    # if ref seq is given, check that it matches some sequence in alignment
        is_valid = False
        for record in msa:
            if refseq in crutil.remove_gaps(record.seq):
                refseq = record.seq
                is_valid = True
                break
        assert is_valid, "reference sequence is not present in given msa"
    else: # if ref seq is not given, use S. cerevisae sequence as ref seq
        for record in msa:
            if record.id == "Saccharomyces":
                refseq = record.seq
                
    
    # adjusting boundary indices
    # **NOTE: this means that all the sites preceding and trailing the last residue
    # in the refseq will never be included, which assumes a high confidence on
    # the quality and reliability of the alignment**
    refseq_gapless = crutil.remove_gaps(refseq)
    left_bound = max(left_bound, 0)
    right_bound = min(right_bound, len(refseq_gapless) - 1)
    
    # converting sequence boundary to alignment boundary
    mapping = crutil.get_mapping(refseq)
    msa_left_bound = mapping[left_bound]
    msa_right_bound = mapping[right_bound]
    
    return msa_left_bound, msa_right_bound


def extract_region_msa(msa, left_bound, right_bound,
                       refseq=None):
    '''
    extract the msa of the target region
    input:
    msa - (alignIO object) msa
    left_bound - (int) left boundary of region in reference sequence
    right_bound - (int) right boundary of region in reference sequence
    refseq - (str, optional) reference sequence
    output:
    region_msa - (alignIO object) msa of target region
    '''
    msa_left_bound, msa_right_bound = get_region_range(msa, left_bound, right_bound,
                                        refseq=refseq)
    region_msa = msa[:, msa_left_bound:(msa_right_bound+1)]
    return region_msa


def extract_random_region_from_msa(msa, regionlen):
    '''
    extract a randomly selected region of specified length from msa
    input:
    msa - (alignIO object) msa
    regionlen - (int) desired region length
    output:
    random_msa - (alignIO object) randomly-selected region msa
    '''
    if len(msa[0].seq) < regionlen:
        return msa
    else:
        starti = random.randrange(len(msa[0].seq) - regionlen)
        random_msa = msa[:, starti:(starti + regionlen)]
        return random_msa


'''
main functions - compute "gappiness" as a measure of alignment quality
algorithm:
- delete sequences with only gaps
- delete longest sequence
- delete shortest sequence
- realignment (remove columns with only gaps)
- compute average number of gaps normalized by length
'''

def get_seq_len(msa):
    '''
    helper function for remove_min_max.
    input:
    msa - (AlignIO object)
    output:
    seqlen - (dict) mapping sequence description to sequence length
    '''
    seqlen = {}
    for record in msa:
        seq_gapless = crutil.remove_gaps(record.seq)
        l = len(seq_gapless)
        desc = record.description
        seqlen[desc] = l
    return seqlen


def filter_msa(msa):
    '''
    Removes the longest sequence, the shortest sequene, and the sequences
    with only gaps from the msa
    input:
    msa - (alignIO object) msa
    output:
    filter_msa - (list of SeqRecord objects) filtered msa
    **NOTE: the output is not an alignIO object**
    '''
    seqlen = get_seq_len(msa)
    seqlen_pos = {desc:l for desc,l in seqlen.items() if l != 0} # excluding 0
    maxlen = max(seqlen_pos, key=seqlen.get)
    minlen = min(seqlen_pos, key=seqlen.get)
    filtered_msa = []
    for record in msa:
        if (record.description != maxlen) and (record.description != minlen)\
            and (seqlen[record.description] != 0):
            filtered_msa.append(record)
        #else:
            #print(record.description)
    return filtered_msa


def get_gap_only_indices(filtered_msa):
    '''
    find the positions that only contain white gaps.
    input:
    filtered_msa - (list of SeqRecord objects)
    output:
    gap_only_indices - (list of int) list of indices corresponding to columns
    that only contain white gaps
    '''
    gap_only_indices = []
    for i in range(len(filtered_msa[0].seq)):
        all_gaps = True
        for record in filtered_msa:
            if record.seq[i] != '-':
                all_gaps = False
                break
        if all_gaps:
            gap_only_indices.append(i)
    return gap_only_indices


def remove_gap_only_indices(filtered_msa):
    '''
    remove the columns that only contain white gaps.
    **NOTE: this function modifies SeqRecord IN PLACE. use with caution.**
    input:
    filtered_msa - (list of SeqRecord objects)
    '''
    goi = get_gap_only_indices(filtered_msa)
    n = len(filtered_msa[0].seq)
    for record in filtered_msa:
        seq = ''
        for i in range(n):
            if i not in goi:
                seq += record.seq[i]
        record.seq = seq


def get_gap_frequency_per_seq(filtered_msa):
    '''
    compute the average frequency of gaps among all sequences in msa
    input:
    filtered_msa - (list of SeqRecord objects)
    output:
    gapfreq - (float) average gap frequency
    '''
    gapcnt = []
    for record in filtered_msa:
        cnt = 0
        for aa in record.seq:
            if aa == '-':
                cnt += 1
        gapcnt.append(cnt)
    gapfreq = np.mean(gapcnt) / len(filtered_msa[0].seq)
    return gapfreq


def compute_alignment_quality(region_msa):
    '''
    finally, function for computing alignment quality.
    input:
    region_msa - (alignIO object) msa of target region
    output:
    gapfreq - (float) average gap frequency
    '''
    filtered_msa = filter_msa(region_msa)
    if len(filtered_msa) >= 1:
        remove_gap_only_indices(filtered_msa)
        gapfreq = get_gap_frequency_per_seq(filtered_msa)
    else:
        gapfreq = None
    return gapfreq


'''
functions for computing sequence heterogeneity
'''

def count_aas_in_column(msa, position, states):
    '''
    count aas in a given column in msa
    input:
    msa - (alignIO object)
    position - (int) the index of the column in msa
    states - (list) all possible aas
    output:
    aa_count - (dict) raw count for each aa in states
    '''
    aa_count = {}
    for state in states:
        aa_count[state] = 0
    for seq in msa:
        aa = seq.seq[position]
        if aa in states:
            aa_count[aa] += 1
    return aa_count


def get_column_aa_distribution(msa, position, states):
    '''
    compute the first order probabilistic distribution of aas in given column
    input:
    msa - (alignIO object)
    position - (int) the index of the column in msa
    states - (list) all possible aas
    output:
    aa_dist - (list) first order probabilistic distribution of aas
    '''
    aa_cnt = count_aas_in_column(msa, position, states)
    total_cnt = sum(aa_cnt.values())
    if total_cnt != 0:
        aa_dist = []
        for k,v in aa_cnt.items():
            aa_dist.append(v / total_cnt)
    else:
        aa_dist = None
    return aa_dist


def get_seq_divergence(msa, states):
    '''
    compute sequence entropy as a measure of heterogenity
    input:
    msa - (alignIO object)
    states - (list) all possible aas
    output:
    S - (float) mean entropy across all columns in msa
    '''
    entropy_list = []
    for p in range(len(msa[0].seq)):
        aa_dist = get_column_aa_distribution(msa, p, states)
        if not aa_dist is None:
            entropy_list.append(stats.entropy(aa_dist))
    S = sum(entropy_list) / len(entropy_list)
    return S


def get_seq_divergence_from_region_msa(region_msa, states, min_n=1):
    '''
    compute sequence divergence for the msa of target region
    input:
    region_msa - (alignIO object) msa of target region
    states - (list) all possible aas
    output:
    seq_div - (float) sequence divergence
    '''
    filtered_msa = filter_msa(region_msa)
    if len(filtered_msa) >= min_n:
        remove_gap_only_indices(filtered_msa)
        seq_div = get_seq_divergence(filtered_msa, states)
    else:
        seq_div = None
    return seq_div