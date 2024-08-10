#!/usr/bin/env python
# coding: utf-8
# Author: Brandon D. Jordan
# Date started: 01/2024
# Obj: print part of analyses of blastn query searching for soybean centromeric repeats

# import modules
import pandas as pd
import sys, getopt

def main(argv):
    print(sys.argv[1])

    # load blast results in tabulated format; skip header rows.
    gen = pd.read_csv(sys.argv[1], header=None,names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                  'qend', 'sstart', 'send', 'evalue', 'bitscore'],
                  sep='\t', comment='#')

    # sort df by sseqid starts
    gen = gen.sort_values(by = ['sseqid','sstart'])

    # Create list of values [Gm01,Gm02,Gm03...Gm20] or [Chr01...Chr20]

    # create a list of integers 1-9
    chr_list1 = [*range(1,10)] # note what occurs when run without the '*'
    chr_list1 = [*range(1,10)] # create a list of integers 1-9
    chr_list = [str(x) for x in chr_list1] # make the list elements strings
    chr_list = ["0" + s for s in chr_list] # append a '0' to the beginning of each element
    chr_list2 = [*range(10,21)] # create a list of integers 10-20
    chr_list = chr_list+chr_list2 # concatenate the two lists
    chr_list = [str(x) for x in chr_list] # make the list elements strings
    chr_list = ["Chr" + s for s in chr_list] # append 'Gm' to each element of list # append 'Chr' for other cases...

    # subset only centromeric elements found on chr's not scaffolds
    gen_sub = gen[gen.sseqid.str.contains("Chr", regex= True, na=False)].sort_values(by = ['sseqid','sstart'])

    print("Whole genome value cts for "+sys.argv[1])
    print(gen_sub.qseqid.value_counts())

    for chr in chr_list:
        sub2 = gen_sub[(gen_sub.pident>=90)& (gen_sub.length.between(65,174))&
               (gen_sub.sseqid.str.contains(chr, regex= True, na=False))]
        sub2.to_csv('out.'+chr+'.tsv', sep="\t",header=False,index=False)
        print("Subset: Chromosomal value cts for "+sys.argv[1])
        print(sub2.qseqid.value_counts())

    del gen, gen_sub, sub2

if __name__ == "__main__":
    main(sys.argv[1:])



