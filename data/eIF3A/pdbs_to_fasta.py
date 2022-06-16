#! python

from gettext import find
import glob
import gzip
import argparse
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


"""
Run this script to create one fasta-formatted file of the sequences of all the pdb files in a directory.
By default the filename/path will be the name in the resulting PDB file
An option can be specified to look for a species name and make that the name instead.
If this option is selected it is assumed that the file is formatted as an Alphafold2 output PDB: use with caution.


"""

def is_gz_file(filepath):
	"""
	This function was sourced from Stack Overflow:
	https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
	"""
	with open(filepath, 'rb') as test_f:
		return test_f.read(2) == b'\x1f\x8b'


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
                if line[0:8] == b"SEQRES  ":
                    resids+=line[19:71].decode("utf-8").strip().split()
        
        return(''.join([three_to_one_map[aa] for aa in resids]))
    
    else:
        with open(filepath) as file:
            resids = []
            pLDDTs = []
            for line in file:
                if line[0:8] == "SEQRES  ":
                    resids+=line[19:71].strip().split()
        
        return(''.join([three_to_one_map[aa] for aa in resids]))


def find_species_name(filepath):
    """
    Parse a line of a pdb file and look for a species name in the header
    """
    search_string = re.compile("ORGANISM_SCIENTIFIC: (.*)")
    if is_gz_file(filepath):
        with gzip.open(filepath) as file:
            resids = []
            pLDDTs = []
            for line in file:
                if line[0:8] == b"SOURCE  ":
                    m = re.search(search_string, line)
                    if m:
                        cleaned = re.sub(r'[^A-Za-z0-9 ]+', '', m.group(1)).strip()
                        # replace whitespace with underscore; muscle breaks at whitespace leading to duplicate entries
                        cleaned = re.sub(" ", "_", cleaned)
                        return(cleaned)
    
    else:
        with open(filepath) as file:
            resids = []
            pLDDTs = []
            for line in file:
                if line[0:8] == "SOURCE  ":
                    m = re.search(search_string, line)
                    if m:
                        cleaned = re.sub(r'[^A-Za-z0-9 ]+', '', m.group(1)).strip()
                        # replace whitespace with underscore; muscle breaks at whitespace leading to duplicate entries
                        cleaned = re.sub(" ", "_", cleaned)
                        return(cleaned)



if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Save sequences to fasta")
	# Required arguments
    parser.add_argument(dest="fpath", help="path to directory with pdb files")
    # Optional arguments
    parser.add_argument("--out_dir", dest="out_directory", default=None, help="Filepath to output directory")
    parser.add_argument("--out", dest="out_fname", default=None, help="output FASTA filename")
    parser.add_argument("--rename", action='store_true', help="Add this flag to attempt to rename species in output fasta")
    options = parser.parse_args()
    
    files = glob.glob(options.fpath + "/*.pdb")
    
    seqs = []
    
    for pdb in files:
        pdb_seq = Seq(read_seq_from_pdb(pdb))
        if options.rename:
            species_name = find_species_name(pdb)
            print(species_name)
        else:
            species_name = pdb
        rec = SeqRecord(pdb_seq, id = species_name, description=pdb)
        seqs.append(rec)

	#print(seqs)
    if not options.out_directory:
        if not options.out_fname:
            SeqIO.write(seqs, "./seqs.fa", "fasta")
        else:
            SeqIO.write(seqs, "./{}".format(options.out_fname), "fasta")
    else:
        if not options.out_fname:
            SeqIO.write(seqs, options.out_directory+"seqs.fa", "fasta")
        else:
            SeqIO.write(seqs, options.out_directory+options.out_fname, "fasta")