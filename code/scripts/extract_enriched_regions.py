import sys, argparse, csv
from Bio import SeqIO
from datetime import datetime
import fractional_charge as fc
import charge_distribution as cd

'''
to write in:
- systematic name/id
- description
- region boundaries
- region sequence
- region length
- fraction in first half of alphabet
'''

if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Extract regions from yeast transcriptome that have a high frequency of letters in the first half of the alphabet")
    parser.add_argument("transcriptome", help="fasta file containing the transcriptome")
    parser.add_argument("target_aas", type=str, help="Look for enrichment in these amino acids. Format = 'AED' for enrichment of Alanine, Glutamate and Aspartate.")
    parser.add_argument("win_size", type=int, help="window size for sliding window analysis")
    parser.add_argument("threshold", type=float, help="cutoff for fraction of sequence enriched in target amino acids")
    parser.add_argument("min_len", type=int, help="minimum region length")
    parser.add_argument("tol", type=int, help="region detection tolerance")
    parser.add_argument("-o", "--output", help="output filename")
    args = parser.parse_args()

    target_list = [i for i in args.target_aas]
    print("Detecting regions enriched in amino acids {}".format(args.target_aas))

    orfs = []
    for record in SeqIO.parse(args.transcriptome, "fasta"):
        orfs.append(record)

    outFile = open(args.output, 'w')
    with outFile:
        outFile.write("# Run started {}\n".format(datetime.now()))
        outFile.write("# Command: {}\n".format(' '.join(sys.argv)))
        outFile.write("# Parameters:\n")
        argsdict = vars(args)
        for k, v in argsdict.items():
            outFile.write("#\t{k}: {v}\n".format(k=k, v=v))

        outHeaders = ['orf', 'gene', 'seq.len', 'left.bound', 'right.bound', 'region.seq', 'region.len',
                      'fraction.enriched.{}'.format(args.target_aas)]
        writer = csv.DictWriter(outFile, fieldnames=outHeaders)
        writer.writeheader()
        for orf in orfs:
            region_boundaries = fc.find_enriched_start_end(orf.seq, target_list, args.win_size,\
                                    args.threshold, args.min_len, args.tol)
            for i in region_boundaries:
                rv = {}
                rv['orf'] = orf.name

                desc = ''
                for c in orf.description:
                    if c != ',':
                        desc += c
                    else:
                        break
                gene = desc.split()[1]
                rv['gene'] = gene

                region_seq = fc.range_to_seq(orf.seq, i)
                rv['seq.len'] = len(orf.seq)
                rv['left.bound'] = i[0]
                rv['right.bound'] = i[1]
                rv['region.seq'] = region_seq
                rv['region.len'] = len(region_seq)
                rv['fraction.enriched.{}'.format(args.target_aas)] = fc.get_enriched_fraction(region_seq, target_list)
                writer.writerow(rv)
        outFile.write("# Run finished {}\n".format(datetime.now()))