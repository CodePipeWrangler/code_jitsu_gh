#!/usr/bin/python3
#Author: Brandon D. Jordan
#Obj: add sequence designation of text and incrementing digit to every title line of fasta file.

import sys, getopt

def main(argv):
   inputfile = ''
   outputfile = ''
   identity = ''

   try:
       opts, args = getopt.getopt(argv,"hi:o:f:",["ifile=", "ofile=", "id="])
   except getopt.GetoptError:
      print ('file.py -i <inputfile> -o <outputfile> -f <identifier>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('file.py -i <inputfile> -o <outputfile> -f <identifier>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
         print ('Output:', outputfile)
      elif opt in ("-f", "--id"):
         identity = arg

   fin = open(inputfile, 'r')
   fout = open(outputfile, 'w')

   linesin = fin.read().splitlines()

   seqnum = 1
   current = None
   for lin in linesin:
     current = lin
     if lin == '>':
        current += identity + '.rpt.' + str(seqnum)
        seqnum = seqnum + 1
     fout.write(current+'\n')

   fout.close()
   fin.close()

if __name__ == "__main__":
    main(sys.argv[1:])
