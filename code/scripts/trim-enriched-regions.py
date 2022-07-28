from os import scandir
import sys, argparse, csv
import glob
import pandas as pd
import fractional_charge as fc


def trim_region(seq, target_list):
    n_lefttrim, n_righttrim = 0, 0
    for AA in seq:
        if AA not in target_list:
            n_lefttrim += 1
        else:
            break
    for AA in seq[::-1]:
        if AA not in target_list:
            n_righttrim += 1
        else:
            break
    trimmed_seq = seq[n_lefttrim:(len(seq) - n_righttrim)]
    return trimmed_seq, n_lefttrim, n_righttrim

if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Extract regions from yeast transcriptome that have a high frequency of letters in the first half of the alphabet")
    parser.add_argument("input_dir", help="path to data; data should be formatted as csvs with the last column being the enriched fraction and the last four characters of the column name are the four enriched amino acids")
    parser.add_argument("output_dir", help="path to output trimmed regions")
    args = parser.parse_args()

    files = glob.glob(args.input_dir+"*.csv")
    
    for file in files:
        df = pd.read_csv(file, comment="#")

        colnames = list(df.columns)
        targets = list(colnames[-1][-4:])
        
        dict_trimmed = []
        for index, row in df.iterrows():
            rv = {}
            for name in colnames:
                rv[name] = row[name]
            oldseq = rv['region.seq']
            newseq, n_lefttrim, n_righttrim = trim_region(oldseq, targets)
            rv['left.bound'] = row['left.bound'] + n_lefttrim
            rv['right.bound'] = row['right.bound'] - n_righttrim
            rv['region.seq'] = newseq
            rv['region.len'] = row['region.len'] - (n_lefttrim + n_righttrim)
            rv['fraction.enriched.{}'.format(''.join(targets))] = fc.get_enriched_fraction(newseq, targets)
            dict_trimmed.append(rv)
        df_trimmed = pd.DataFrame.from_records(dict_trimmed)

        # print(df.tail())
        # print(df_trimmed.tail())
        # print(oldseq)
        # print(newseq)

        df_trimmed.to_csv(args.output_dir+'/'+''.join(targets)+"_enriched_trimmed.csv")