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
- fractional charge
- (kappa)
'''

if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Extract charged regions from yeast transcriptome")
    parser.add_argument("transcriptome", help="fasta file containing the transcriptome")
    parser.add_argument("win_size", type=int, help="window size for sliding window analysis")
    parser.add_argument("threshold", type=float, help="charge threshold")
    parser.add_argument("min_len", type=int, help="minimum region length")
    parser.add_argument("tol", type=int, help="region detection tolerance")
    parser.add_argument("-o", "--output", help="output filename")
    args = parser.parse_args()

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
                      'charge.asymmetry', 'frac.charge', 'kappa1', 'kappa2']
        writer = csv.DictWriter(outFile, fieldnames=outHeaders)
        writer.writeheader()
        for orf in orfs:
            region_boundaries = fc.find_start_end(orf.seq, args.win_size,\
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
                rv['charge.asymmetry'] = cd.get_sigma(region_seq)
                rv['frac.charge'] = fc.get_fractional_charge(region_seq)
                rv['kappa1'] = cd.get_avg_kappa(region_seq, alt_norm=False)
                rv['kappa2'] = cd.get_avg_kappa(region_seq, alt_norm=True)
                writer.writerow(rv)
        outFile.write("# Run finished {}\n".format(datetime.now()))