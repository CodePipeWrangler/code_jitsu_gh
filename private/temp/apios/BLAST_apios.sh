# BLAST WORKFLOW APIOS

cd WORKING_DIR
mkdir blastdb blastout query

APIDATA=/pscratch/sd/j/jordan25/genomics/apios
GNM_apiam=$APIDATA/apiam.LA2127.gnm1.KVMK/apiam.LA2127.gnm1.KVMK.genome_main.fna
GNM_apipr=$APIDATA/apipr.MO19963523.gnm1.814V/apipr.MO19963523.gnm1.814V.genome_main.fna

cat $GNM_apiam | makeblastdb -in - -dbtype nucl -title Apiam -hash_index -parse_seqids -out blastdb/Apiam &
cat $GNM_apipr | makeblastdb -in - -dbtype nucl -title Apipr -hash_index -parse_seqids -out blastdb/Apipr &

I went back and made more consensus sequences to query from the alignments of 193 bp repeat sequences and also from apiam. Possibly sone overalap, but there are 9 consensus sequences total.
Seven of these are from apipr, and two are from apiam.

# Starting query
>Type_1.1.apipr.cons
AAAAACATTTAAAATATTATGCnACGAGAAAGTACAGTTTTATAAAAATGTAAAACCGTT
TCAAGAGCTTAAAGTTTCAGAATTTCAGCATACTAAATCTGTGTTACTAAATATATTAAA
AATAGTAAAAGAAGTCATAAAATAACTTATTATATGTCAAATTAAAGCTTAGTATGTCTA
TTTTATGTTTATGTn
>Type_1.2.apipr.cons
nnAAAACTGTACTATTTTGGTTATAATAAGTTAAATGTTTTTACATAACAATTAAATAGA
CATCTTAAGCTTTGATTTGACACATTATAAGCTATATTATGACTTCGTTTACTATTTTTA
ATATATTTCTTAAGACAGATTTAATGTATTAAAATTCTGAAAATTTATATGCTGAAACAA
TTTTAGTATTTTTGG
>Type_1.3.apipr.cons
AAAAATAAAATAGACATCTTGAGCTTTAATTTGACATATAATAATCTATTTTATGATTTC
GTTTACTAATTTTAATATATTTATCAAAACAGATTTATTGTATTAGAATTCTGAAAATAT
CTGTTGTGCAACTGTTTTTATAAATTTAGAAAACTGTACTACCTCAGTTATAATAATTTA
AATGTTTTTGTCT
>Type_2.apipr.cons
AAAAATAGTAAACAAAGTCATAAAATAGCTTATAATATGTCAAATAAAAGCTTATGATGT
CTATTTTATTTTTATGTAAAAATATTTAACnTATTATAAnGGAAATAGTACAGTTTTACA
AAATTCATAAAATTGTTCCAGAATCTGAATTTTCAGAATTTCAGTACGGTAAATTTGGTT
TAATAAATATATT
>Type_3.apipr.cons
AAAAAAAGTCATAAAATAGCTTATAATATGTCAAATAAAAGCTTATGATGTCTATTTTAT
TGTTATGTAAAAATATTTAACGTATTATAACGGAAATAGTACAGTTTTACAAAATTCATA
AAATTGTTCCAGAATCTGAATTTTCAGAATTTCAGTACGGTAAATTTGGTTTAATAAATA
TATTAAAAATAGT
>Type_4.apipr.cons
AAAACAGATTTATTGTATTAGAATTCTGAAAATTTCTGTTGTGCAACTGATTTTATAATT
TTGGAAAACTnnnGTACTACCTTGTTTATAATAATTTAAATATTTTCACCTAATAATGAA
AnTAGACATCTTAAGCTTTAATTTGACATATTACAAGTTATTTTATGAATTCGTTTACCA
TTTTTATTATATTTCTTnnnnnnnnnnnnnnnnn
>Type_5.apipr.cons
AAAACATTTTaTAAATTTTTGAAAACTTTTCTTTTTTaCTTTTTATTACTTTAATTTTTT
TTCCTTACCATTAAATTAACCTTCTTAACTTTTTTTTTACCATTTTTAAACTTTTTTTTT
ACTTTTTTTTCTTTTTTTTATTTTTTTTTTTAAACCAATTTTATCTTATTAAAATTTTAA
AATTTAAATTTTT
>Type_1.apiam.cons
AAAACTTTTCTTTTTTGGTTATTATTACTTTAATTTTTTTTCCTTACCATTAAATTAACC
TTCTTAACTTTTATTTTACCTTTTTTTAACTTTTTTTTTACTTTGTTTTCTTTTTTTTAT
TTTTTTTTTTAAACCAATTTTACCTTATTAAATTTTTAAnAATTTAAGGTTTTGAAAACT
TTTTTTTATTTTTG
>Type_2.apiam.cons
nAAAAAATTTTACCTTTTTTTAACCAAATTATTCCATTTTTCCAAAATCATAAACATTTT
TCCAAACCTTAATTTTTCAAAATTTTATTACGTTAATTTTGTTTTAAAAATTTTTTTAAA
ATTATTAACCAAATTCTTAAATTAGTTTTTATTTTTTCAATTCAAAGTTTAGATTTTTTT
TTTATTTTTTTTTT


# BLAST searches
  for path in blastdb/*.nsq; do 
    file=`basename $path .nsq`
    echo "WORKING on $file"
    blastn -query query/apipr_ALL.cons -db blastdb/$file -out blastout/cent_rpt.x.$file.bln \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" &
  done
  
# Find the most frequent target sequences
  cat blastout/cent_rpt.x.Apiam.bln | awk '$3>=90 {print $13}' | sort | uniq -c | sort -n | tail -1 | awk '{print $2}'
  # Removing the length ($4) constraint was necessary to obtain results because matches were no greater than 69 bp
    >Apiam_freq
    AAAACTGTACTATCTTGGTTATAATAAGTTAAATGTTTTTACATAACAATTAAATAGACATCTTAAGCTTTGATTTGACATATTATAAGCTATTTTATGACTTCGTTTACTATTTTTAATATATTTCTTAAGACAGATTTAATGTATTAAAATTCTGAAAATTTATATTCTGAAACAGTTTTAGTA

  cat blastout/cent_rpt.x.Apipr.bln | awk '$3>=90 && $4==193 {print $13}' | sort | uniq -c | sort -n | tail -1 | awk '{print $2}'
  # Removing the length ($4) constraint was necessary to obtain results because matches were no greater than 69 bp
    >Apipr_freq
    AAAAATAGTAAACGAAGTCATAATATAGCTTATAATGTGTCAAATCAAAGCTTAAGATGTCTATTTAATTGTTATGTAAAAACATTTAACTTATTATAACCAAAATAGTACAGTTTTCCAAAAT-CATAAAATTGTTCCAGAATCTGAATTTTCAGAATTTCAGTACGGTAAATTTGGTTTAATAAATATATT-CATAAAATTGTTCCAGAATCTGAATTTTCAGAATTTCAGTACGGTAAATTTGGTTTAATAAATATATT


# Here is the distance matrix for the consensus sequences created via clustalo
Type_1.1.apipr.cons 100.000000 29.230769 30.569948 38.860104 43.523316 33.846154 25.906736 29.896907 48.969072 30.107527 37.823834
Type_1.2.apipr.cons 29.230769 100.000000 60.103627 32.124352 33.678756 52.307692 50.259067 56.185567 36.082474 97.311828 32.642487
Type_1.3.apipr.cons 30.569948 60.103627 100.000000 33.160622 34.196891 47.668394 34.715026 35.751295 30.051813 63.440860 32.124352
Type_2.apipr.cons   38.860104 32.124352 33.160622 100.000000 93.264249 30.569948 30.569948 28.497409 36.787565 34.408602 92.746114
Type_3.apipr.cons   43.523316 33.678756 34.196891 93.264249 100.000000 31.606218 33.678756 30.051813 40.414508 34.408602 88.601036
Type_4.apipr.cons   33.846154 52.307692 47.668394 30.569948 31.606218 100.000000 43.005181 37.113402 35.051546 56.451613 29.015544
Type_5.apipr.cons   25.906736 50.259067 34.715026 30.569948 33.678756 43.005181 100.000000 81.347150 41.968912 52.688172 29.015544
Type_1.apiam.cons   29.896907 56.185567 35.751295 28.497409 30.051813 37.113402 81.347150 100.000000 43.814433 56.451613 28.497409
Type_2.apiam.cons   48.969072 36.082474 30.051813 36.787565 40.414508 35.051546 41.968912 43.814433 100.000000 37.096774 40.414508
Apiam_freq          30.107527 97.311828 63.440860 34.408602 34.408602 56.451613 52.688172 56.451613 37.096774 100.000000 33.870968
Apipr_freq          37.823834 32.642487 32.124352 92.746114 88.601036 29.015544 29.015544 28.497409 40.414508 33.870968 100.000000

# From here I removed consensus sequences that were similar to one another (>=80% identity) because I can adjust the filtering setting s to capture these during
# data analysis. Type_1.2.apipr.cons was >97% identical to Apiam_freq, so I opted for the latter as the representative sequence of the two. Type_2.apipr.cons and 
# Type_3.apipr.cons were similar to Apipr_freq, so I opted for the latter as the representative sequence of the three. Type_5.apipr.cons and Type_1.apiam.cons were
# similar but from different species so i kept both. 

# This leaves the following consensus sequences to query against the blast databases:
cat query/apios_ALL.cons
>Type_1.1.apipr.cons
AAAAACATTTAAAATATTATGCnACGAGAAAGTACAGTTTTATAAAAATGTAAAACCGTT
TCAAGAGCTTAAAGTTTCAGAATTTCAGCATACTAAATCTGTGTTACTAAATATATTAAA
AATAGTAAAAGAAGTCATAAAATAACTTATTATATGTCAAATTAAAGCTTAGTATGTCTA
TTTTATGTTTATGTn
>Apiam_freq
AAAACTGTACTATCTTGGTTATAATAAGTTAAATGTTTTTACATAACAATTAAATAGACATCTTAAGCTTTGATTTGACATATTATAAGCTATTTTATGACTTCGTTTACTATTTTTAATATATTTCTTAAGACAGATTTAATGTATTAAAATTCTGAAAATTTATATTCTGAAACAGTTTTAGTA
>Type_1.3.apipr.cons
AAAAATAAAATAGACATCTTGAGCTTTAATTTGACATATAATAATCTATTTTATGATTTC
GTTTACTAATTTTAATATATTTATCAAAACAGATTTATTGTATTAGAATTCTGAAAATAT
CTGTTGTGCAACTGTTTTTATAAATTTAGAAAACTGTACTACCTCAGTTATAATAATTTA
AATGTTTTTGTCT
>Apipr_freq
AAAAATAGTAAACGAAGTCATAATATAGCTTATAATGTGTCAAATCAAAGCTTAAGATGTCTATTTAATTGTTATGTAAAAACATTTAACTTATTATAACCAAAATAGTACAGTTTTCCAAAAT-CATAAAATTGTTCCAGAATCTGAATTTTCAGAATTTCAGTACGGTAAATTTGGTTTAATAAATATATT-CATAAAATTGTTCCAGAATCTGAATTTTCAGAATTTCAGTACGGTAAATTTGGTTTAATAAATATATT
>Type_4.apipr.cons
AAAACAGATTTATTGTATTAGAATTCTGAAAATTTCTGTTGTGCAACTGATTTTATAATT
TTGGAAAACTnnnGTACTACCTTGTTTATAATAATTTAAATATTTTCACCTAATAATGAA
AnTAGACATCTTAAGCTTTAATTTGACATATTACAAGTTATTTTATGAATTCGTTTACCA
TTTTTATTATATTTCTTnnnnnnnnnnnnnnnnn
>Type_5.apipr.cons
AAAACATTTTaTAAATTTTTGAAAACTTTTCTTTTTTaCTTTTTATTACTTTAATTTTTT
TTCCTTACCATTAAATTAACCTTCTTAACTTTTTTTTTACCATTTTTAAACTTTTTTTTT
ACTTTTTTTTCTTTTTTTTATTTTTTTTTTTAAACCAATTTTATCTTATTAAAATTTTAA
AATTTAAATTTTT
>Type_1.apiam.cons
AAAACTTTTCTTTTTTGGTTATTATTACTTTAATTTTTTTTCCTTACCATTAAATTAACC
TTCTTAACTTTTATTTTACCTTTTTTTAACTTTTTTTTTACTTTGTTTTCTTTTTTTTAT
TTTTTTTTTTAAACCAATTTTACCTTATTAAATTTTTAAnAATTTAAGGTTTTGAAAACT
TTTTTTTATTTTTG
>Type_2.apiam.cons
nAAAAAATTTTACCTTTTTTTAACCAAATTATTCCATTTTTCCAAAATCATAAACATTTT
TCCAAACCTTAATTTTTCAAAATTTTATTACGTTAATTTTGTTTTAAAAATTTTTTTAAA
ATTATTAACCAAATTCTTAAATTAGTTTTTATTTTTTCAATTCAAAGTTTAGATTTTTTT
TTTATTTTTTTTTT

# BLAST searches
  for path in blastdb/*.nsq; do 
    file=`basename $path .nsq`
    echo "WORKING on $file"
    blastn -query query/apipr_ALL.cons -db blastdb/$file -out blastout/cent_rpt.x.$file.bln \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" &
  done

# Summarize results
  for filepath in blastout/*bln; do 
    echo $filepath
    cut -f1 $filepath | sort | uniq -c;
    echo
  done

<< RESULT
  blastout/cent_rpt.x.Apaim.bln
   1040 Apiam_freq
    296 Apipr_freq
    499 Type_1.1.apipr.cons
    667 Type_1.3.apipr.cons
    236 Type_4.apipr.cons

  blastout/cent_rpt.x.Apipr.bln
   2279 Apiam_freq
    761 Apipr_freq
   3096 Type_1.1.apipr.cons
   3896 Type_1.3.apipr.cons
   1530 Type_4.apipr.cons

# Extract coords into bed files
  cat blastout/cent_rpt.x.Apiam.bln | 
    awk -v OFS="\t" '$10>$9 {print $2, $9-1, $10}
                     $10<$9 {print $2, $10, $9-1}' | 
    sort -k1,1 -k2n,2n > blastout/cent_rpt.x.Apiam.bed

  cat blastout/cent_rpt.x.Apipr.bln |
    awk -v OFS="\t" '$10>$9 {print $2, $9-1, $10}
                     $10<$9 {print $2, $10, $9-1}' | 
    sort -k1,1 -k2n,2n > blastout/cent_rpt.x.Apipr.bed

# Merge hits within specified distance of one another
  for width in 10 25 50 100; do
    echo "Merge width: ${width}000"
    bedtools merge -i blastout/cent_rpt.x.Apiam.bed -d ${width}000 > blastout/Apiam_merged_${width}k.bed
    bedtools merge -i blastout/cent_rpt.x.Apipr.bed -d ${width}000 > blastout/Apipr_merged_${width}k.bed
  done

  for filepath in blastout/*merged*; do 
    base=`basename $filepath`
    echo $base
    awk -v OFS="\t" 'NR==1 {span=$3-$2; prev=$1; grand+=$3-$2} 
                     NR>1 && $1==prev {span+=$3-$2; grand+=$3-$2} 
                     NR>1 && $1!=prev {print prev, span; prev=$1; span=$3-$2; grand+=$3-$2} 
                     END{print $1, span; print "Total:", grand}' $filepath
    echo
  done

<< RESULT

Apiam_merged_100k.bed
apiam.LA2127.gnm1.Chr01 834710
apiam.LA2127.gnm1.Chr02 1753293
apiam.LA2127.gnm1.Chr03 208661
apiam.LA2127.gnm1.Chr04 1950844
apiam.LA2127.gnm1.Chr05 774025
apiam.LA2127.gnm1.Chr06 276292
apiam.LA2127.gnm1.Chr07 2069123
apiam.LA2127.gnm1.Chr08 609417
apiam.LA2127.gnm1.Chr09 4107603
apiam.LA2127.gnm1.Chr10 1295598
apiam.LA2127.gnm1.Chr11 22473
Total:  13902039

Apiam_merged_10k.bed
apiam.LA2127.gnm1.Chr01 233829
apiam.LA2127.gnm1.Chr02 269143
apiam.LA2127.gnm1.Chr03 83178
apiam.LA2127.gnm1.Chr04 261265
apiam.LA2127.gnm1.Chr05 106332
apiam.LA2127.gnm1.Chr06 105532
apiam.LA2127.gnm1.Chr07 281902
apiam.LA2127.gnm1.Chr08 117258
apiam.LA2127.gnm1.Chr09 334126
apiam.LA2127.gnm1.Chr10 143666
apiam.LA2127.gnm1.Chr11 22473
Total:  1958704

Apiam_merged_25k.bed
apiam.LA2127.gnm1.Chr01 278599
apiam.LA2127.gnm1.Chr02 464313
apiam.LA2127.gnm1.Chr03 144238
apiam.LA2127.gnm1.Chr04 629573
apiam.LA2127.gnm1.Chr05 203990
apiam.LA2127.gnm1.Chr06 151500
apiam.LA2127.gnm1.Chr07 681204
apiam.LA2127.gnm1.Chr08 185183
apiam.LA2127.gnm1.Chr09 645542
apiam.LA2127.gnm1.Chr10 330759
apiam.LA2127.gnm1.Chr11 22473
Total:  3737374

Apiam_merged_50k.bed
apiam.LA2127.gnm1.Chr01 690306
apiam.LA2127.gnm1.Chr02 1072164
apiam.LA2127.gnm1.Chr03 144238
apiam.LA2127.gnm1.Chr04 1208069
apiam.LA2127.gnm1.Chr05 241898
apiam.LA2127.gnm1.Chr06 151500
apiam.LA2127.gnm1.Chr07 1208189
apiam.LA2127.gnm1.Chr08 325544
apiam.LA2127.gnm1.Chr09 1274103
apiam.LA2127.gnm1.Chr10 704961
apiam.LA2127.gnm1.Chr11 22473
Total:  7043445

Apipr_merged_100k.bed
apipr.MO19963523.gnm1.Chr01     3331577
apipr.MO19963523.gnm1.Chr02     5086684
apipr.MO19963523.gnm1.Chr03     2845815
apipr.MO19963523.gnm1.Chr04     3265912
apipr.MO19963523.gnm1.Chr05     3068840
apipr.MO19963523.gnm1.Chr06     4006882
apipr.MO19963523.gnm1.Chr07     4112658
apipr.MO19963523.gnm1.Chr08     3482862
apipr.MO19963523.gnm1.Chr09     3082970
apipr.MO19963523.gnm1.Chr10     3501914
apipr.MO19963523.gnm1.Chr11     3501783
Total:  39287897

Apipr_merged_10k.bed
apipr.MO19963523.gnm1.Chr01     948506
apipr.MO19963523.gnm1.Chr02     1095078
apipr.MO19963523.gnm1.Chr03     567922
apipr.MO19963523.gnm1.Chr04     775167
apipr.MO19963523.gnm1.Chr05     685998
apipr.MO19963523.gnm1.Chr06     888324
apipr.MO19963523.gnm1.Chr07     1185506
apipr.MO19963523.gnm1.Chr08     728628
apipr.MO19963523.gnm1.Chr09     581765
apipr.MO19963523.gnm1.Chr10     824637
apipr.MO19963523.gnm1.Chr11     694255
Total:  8975786

Apipr_merged_25k.bed
apipr.MO19963523.gnm1.Chr01     1379396
apipr.MO19963523.gnm1.Chr02     2243135
apipr.MO19963523.gnm1.Chr03     1048856
apipr.MO19963523.gnm1.Chr04     1226640
apipr.MO19963523.gnm1.Chr05     1010539
apipr.MO19963523.gnm1.Chr06     1462299
apipr.MO19963523.gnm1.Chr07     2106443
apipr.MO19963523.gnm1.Chr08     1174309
apipr.MO19963523.gnm1.Chr09     969089
apipr.MO19963523.gnm1.Chr10     1315977
apipr.MO19963523.gnm1.Chr11     1240005
Total:  15176688

Apipr_merged_50k.bed
apipr.MO19963523.gnm1.Chr01     2290055
apipr.MO19963523.gnm1.Chr02     3585396
apipr.MO19963523.gnm1.Chr03     1642809
apipr.MO19963523.gnm1.Chr04     1645633
apipr.MO19963523.gnm1.Chr05     1727310
apipr.MO19963523.gnm1.Chr06     2311071
apipr.MO19963523.gnm1.Chr07     3032142
apipr.MO19963523.gnm1.Chr08     2219172
apipr.MO19963523.gnm1.Chr09     1762056
apipr.MO19963523.gnm1.Chr10     2122790
apipr.MO19963523.gnm1.Chr11     2244818
Total:  24583252

I opted for the merged_50k bed file results
--------------------------------------------------------------------
# BLAST WORKFLOW Gmax

cd WORKING_DIR
mkdir blastdb blastout

GLYMADATA=/pscratch/sd/j/jordan25/genomics/glycine/g.max
GNM_Wm826=$GLYMADATA/glyma.Wm82.gnm6.S97D/glyma.Wm82.gnm6.S97D.genome_main.fna
GNM_Wm825=$GLYMADATA/glyma.Wm82.gnm5.NRKG/glyma.Wm82.gnm5.NRKG.genome_main.fna
GNM_Lee3=$GLYMADATA/glyma.Lee.gnm3.VG1C/glyma.Lee.gnm3.VG1C.genome_main.fna


cat $GNM_Wm826 | makeblastdb -in - -dbtype nucl -title Wm826 -hash_index -parse_seqids -out blastdb/Wm826 &
cat $GNM_Wm825 | makeblastdb -in - -dbtype nucl -title Wm825 -hash_index -parse_seqids -out blastdb/Wm825 &
cat $GNM_Lee3 | makeblastdb -in - -dbtype nucl -title Lee3 -hash_index -parse_seqids -out blastdb/Lee3 &

# Starting query
>CentGm-1
TGTGAAAAGTTATGACCATTTGAATTTCTCGAGAGCTTCCGTTGTTCAATTTCGAGCGTCTCGATATATTATGCGCCTGAATCGGACATCCG
>CentGm-2
AGTCAAAAGTTATTGTCGTTTGACTTTTCTCAGAGCTTCCGTTTTCAATTACGAGCGTCTCGATATATTACGGGACTCAATCGGACATCCG
>Glyma_91bp_cons
AAAAAGTTATTGTCGTTTGAATTTGCTCAGAGCTTCAACATTCAATTTCGAGCGTCTCGATATATTACGGGACTCAATCAGACATCCGAGT
>Glyma_92bp_cons
AAAAGTTATGACCATTCGAATTTCTCGAGAGCTTCCGTTGTTCAATTTCGAGCGTCTCGATATATTATGTCCCCGAATCGGACATTCGTGTG

# BLAST searches
  for path in blastdb/*.nsq; do 
    file=`basename $path .nsq`
    echo "WORKING on $file"
    blastn -query query/glyma.cons -db blastdb/$file -out blastout/cent_rpt.x.$file.bln \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" &
  done

# Summarize results
  for filepath in blastout/*bln; do 
    echo $filepath
    cut -f1 $filepath | sort | uniq -c;
    echo
  done

blastout/cent_rpt.x.Lee3.bln
 149071 CentGm-1
  27755 CentGm-2
  75285 Glyma_91bp_cons
  94480 Glyma_92bp_cons

blastout/cent_rpt.x.Wm825.bln
 127812 CentGm-1
  29529 CentGm-2
  75903 Glyma_91bp_cons
  82147 Glyma_92bp_cons

blastout/cent_rpt.x.Wm826.bln
 127049 CentGm-1
  29295 CentGm-2
  75533 Glyma_91bp_cons

# BLAST WORKFLOW Perr. glycine

cd WORKING_DIR
mkdir blastdb blastout query

PERRGLYDATA=/pscratch/sd/j/jordan25/genomics/glycine/perr_gly
GNM_glycy=$PERRGLYDATA/glycy.G1267.gnm1.YWW6/glycy.G1267.gnm1.YWW6.genome_main.fna
GNM_glyd3=$PERRGLYDATA/glyd3.G1403.gnm1.CL6K/glyd3.G1403.gnm1.CL6K.genome_main.fna
GNM_glydo=$PERRGLYDATA/glydo.G1134.gnm1.PP7B/glydo.G1134.gnm1.PP7B.genome_main.fna
GNM_glyst=$PERRGLYDATA/glyst.G1974.gnm1.7MZB/glyst.G1974.gnm1.7MZB.genome_main.fna
GNM_glysy=$PERRGLYDATA/glysy.G1300.gnm1.C11H/glysy.G1300.gnm1.C11H.genome_main.fna
GNM_glyfa=$PERRGLYDATA/glyfa.G1718.gnm1.B1PY/glyfa.G1718.gnm1.B1PY.genome_main.fna


cat $GNM_glycy | makeblastdb -in - -dbtype nucl -title glycy -hash_index -parse_seqids -out blastdb/glycy &
cat $GNM_glyd3 | makeblastdb -in - -dbtype nucl -title glyd3 -hash_index -parse_seqids -out blastdb/glyd3 &
cat $GNM_glydo | makeblastdb -in - -dbtype nucl -title glydo -hash_index -parse_seqids -out blastdb/glydo &
cat $GNM_glyst | makeblastdb -in - -dbtype nucl -title glyst -hash_index -parse_seqids -out blastdb/glyst &
cat $GNM_glysy | makeblastdb -in - -dbtype nucl -title glysy -hash_index -parse_seqids -out blastdb/glysy &
cat $GNM_glyfa | makeblastdb -in - -dbtype nucl -title glyfa -hash_index -parse_seqids -out blastdb/glyfa &

# Starting query
>CentGm-1
TGTGAAAAGTTATGACCATTTGAATTTCTCGAGAGCTTCCGTTGTTCAATTTCGAGCGTCTCGATATATTATGCGCCTGAATCGGACATCCG
>CentGm-2
AGTCAAAAGTTATTGTCGTTTGACTTTTCTCAGAGCTTCCGTTTTCAATTACGAGCGTCTCGATATATTACGGGACTCAATCGGACATCCG
>Glyma_91bp_cons
AAAAAGTTATTGTCGTTTGAATTTGCTCAGAGCTTCAACATTCAATTTCGAGCGTCTCGATATATTACGGGACTCAATCAGACATCCGAGT
>Glyma_92bp_cons
AAAAGTTATGACCATTCGAATTTCTCGAGAGCTTCCGTTGTTCAATTTCGAGCGTCTCGATATATTATGTCCCCGAATCGGACATTCGTGTG

# BLAST searches
  for path in blastdb/*.nsq; do 
    file=`basename $path .nsq`
    echo "WORKING on $file"
    blastn -query query/glyma.cons -db blastdb/$file -out blastout/cent_rpt.x.$file.bln \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" &
  done

# Summarize results
  for filepath in blastout/*bln; do 
    echo $filepath
    cut -f1 $filepath | sort | uniq -c;
    echo
  done

# blastout/cent_rpt.x.glycy.bln

# blastout/cent_rpt.x.glyd3.bln

# blastout/cent_rpt.x.glydo.bln

# blastout/cent_rpt.x.glyfa.bln
#    5670 CentGm-1
#      16 Glyma_91bp_cons
#     409 Glyma_92bp_cons

# blastout/cent_rpt.x.glyst.bln
#      13 CentGm-1
#       3 CentGm-2
#       4 Glyma_91bp_cons
#      13 Glyma_92bp_cons

# blastout/cent_rpt.x.glysy.bln
#       1 CentGm-1
#       1 Glyma_92bp_cons

# Find the most frequent target sequences
  cat blastout/cent_rpt.x.glyfa.bln | awk '$3>=90 && $4==91 {print $13}' | sort | uniq -c | sort -n | tail -1 | awk '{print $2}'
    >
    

# I went back again and with my new machine learning workflow. I wanted to redo this analysis of the 193 bp repeats.
#

(
for filename in ultra*tsv; do
    ref=`echo $filename | perl -pe 's/ultra\.(.+)\.p\d+\.tsv/$1/'` # EDIT ACCORDING TO YOUR FILENAME # notice the updated code must be integrated into other scripts
    echo -e "\n$ref";
    for x in 193; do
        for i in {01..11}; do
            maxBin=`awk -v CUT="$x" ''/'[Cc]hr'$i'/ && $4==CUT {print int($2/1000000)}' $filename | sort -n | uniq -c | awk -v max=0 '{if($1>max){want=$2; max=$1}}END{print want}'`;
            pos1=$(($maxBin*1000000)); pos2=$((pos1+1000000));
            echo -e "Chr$i\t$x\t$pos1\t$pos2";
            awk -v CUT="$x" -v START="$pos1" -v STOP="$pos2" '$4==CUT && $2>=START && $2<=STOP && '/'[Cc]hr'$i'/ {print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $filename >> $ref.$x.fna ;
        done
    done
done
) &

# mmseqs2 
mkdir clust_mod0 clust_mod1 clust_mod2 clust_mod3 clust_hash.90 db aln prefilter results;
mmseqs createdb *fna db/inputDB
(
for i in {0..3}; do # EDIT THE CLUSTERING MODES OUTPUT HERE
    mmseqs cluster db/inputDB clust_mod$i/clusteredDB tmp --cluster-mode $i --threads 8 && 
    mmseqs createtsv db/inputDB db/inputDB clust_mod$i/clusteredDB clust_mod$i/clustered_pairs.tsv
    mmseqs createsubdb clust_mod$i/clusteredDB db/inputDB db/repsDB
    mmseqs convert2fasta db/repsDB clust_mod$i/reps.fasta
done
)

### find which clusters are the largest (all representative sequences of clusters in rep.fasta can be used to search for homology against the genome of interest)

for file in clust*/clust*tsv; do
    echo "$file --------------------"
    awk '{print $1}' $file | uniq -c | sort -r | head -30
    echo "------------------"
done

# Extract clustered (derived via cluster-mode 1) sequences from the original fasta file using the representative sequences

vim rep_list.txt # contains the IDS for sequences of clusters to be extracted from the fasta file

i=1
while read -r rep; do
  rm -f tmp_ids.txt
  awk -v rep="$rep" '$1==rep {print $2}' clust_mod1/clustered_pairs.tsv > tmp_ids.txt # change target cluster if neeeded
  seqkit grep -f tmp_ids.txt *fna > clust${i}_${rep}_seqs.fasta
  ((i++))
done < rep_list.txt

# cluster grouped filenames and sizes for Apiam
clust_apiam.LA2127.gnm1.Chr01_168551079_586_193_seqs.fasta = clust 1
38
AAAAAATTTTTTTATTTTTGAAAACTTTTCTTTTTTGGTTATTATTACTTTAATTTTTTTTCATTACCATTAAATTAACCTTCTTAACTTTTATTTTACCTTTTTTTAACTTTTTTTTTACTTCTTTTTCTTTTTTTTATTTTTTTTTTTAAACCAATTTTACCTTATTAAATTTTTAAAATTTAAGGTTTTG

clust_apiam.LA2127.gnm1.Chr04_48696269_588_193_seqs.fasta = clust 2
29
AAAAAATTTTACCTTTTTTTAACCAAATTATTCCATTTTTCCAAAATCATAAACATTTTTCCAAACCTTAATTTTTCAAAATTTTATTACGTTAATTTTGTTTTAAAAATTTTTTTAAAATTATTAACCAAATTCTTAAATTAGCTTTTATTTTTTCAATTCAAAGTTTAGATTTTTTTTTTATTTTTTTTTT

clust_apiam.LA2127.gnm1.Chr07_85890353_393_193_seqs.fasta = clust 3
10
AAAAAGATTTAAAATATAATGACCGAGAAAGTACAGTTTTACAAAACCATAAAATTGTTTCCAAAGCCTAAAGTTTTCGAATTTCAGCATACTAAATCTGTATTACAAAATATATTAAAAATAGTAAACAAAGTCATAAAATAACTTATTATATGTCAAATTAAAGATTAGTATGTCTATTTTATCTCTATGT

clust_apiam.LA2127.gnm1.Chr09_35025643_393_193_seqs.fasta = clust 4
13
AAAATAGACATCTTGAGCTTTAATTTGACATATAATATGTTATTTTATGATTTCGTTTACTATTTTTATTATATTTCTTAAGACAGATTTACTACACTATAATTCTGAAAATTTCTGTTGTGTAATCATTTTTATAATTTTGGAAAATTGTACTACCTTGGTTATAATAATTTAAATAGTTTTACCTAATAAT

# cluster grouped filenames and sizes for Apiam
clust6_apipr.MO19963523.gnm1.Chr04_51325955_393_193_seqs.fasta
110
AAAAATAAAATAGACATCTTGAGCTTTAATTTGACATATAATAATCTATTTTATGATTTCGTTTACTAATTTTAATATATTTATCAAAACAGATTTATTGTATTAGAATTCTAAAAATATCTGTTGTGCAACTGTTTTTATAAATTTAGAAAATTATACTACCTCAGTTATAATAATTTAAATGTTTTTGTCT

clust7_apipr.MO19963523.gnm1.Chr01_136008924_394_193_seqs.fasta
74
AAAAACATTTAAAATATCATTCACGAGAAAGTACAGTTTTATAAAAATGTAAACCCGTTTCCAAAGCCTCATGTTTCAGAATTTCAGTATACTAAATCTGTGTTACTAAATATATTAAAAATAGTAAACGAAGTCATAAAATAACTTATTATAAGTCGAATTAAAGCTTAGTATGTCTATTTTATGTTTATGT

clust8_apipr.MO19963523.gnm1.Chr05_92366241_389_193_seqs.fasta
56
AAAACAGATTTAATGTATTAAAATTCTGAAAATTTATATGCTGAAACAATTTTAGTATTTTTGGAAAACTGTACTATTTTGGTTATAATAAGTTAAATGTTTTTACATAACAATTAAATAGACATCTTAAGCTTTGATTTGACACATTATAAGCTATATTATGACTTCGTTTACTATTTTTAATATATTTCTT

clust9_apipr.MO19963523.gnm1.Chr04_51329060_394_193_seqs.fasta
13
AAAACAGATTTATTGTATTAGAATTCTGAAAATTTCTGTTGTGCAACTGATTTTATAATTTTGGAAAACTGTACTACCTTGTTTATAATAATTTAAATATTTTCACCTAATAATGAAATAGACATCTTAAGCTTTAATTTGACATATTACAAGTTATTTTATGAATTCGTTTACCATTTTTATTATATTTCTT

# Convert fasta IDs into tab-separated format
seqkit seq -n clust6_apipr.MO19963523.gnm1.Chr04_51325955_393_193_seqs.fasta | perl -pe 's/^(.*?[Cc]hr\d+)_/\1\t/; s/_/\t/g' > clust_6.tsv
seqkit seq -n clust7_apipr.MO19963523.gnm1.Chr01_136008924_394_193_seqs.fasta | perl -pe 's/^(.*?[Cc]hr\d+)_/\1\t/; s/_/\t/g' > clust_7.tsv
seqkit seq -n clust8_apipr.MO19963523.gnm1.Chr05_92366241_389_193_seqs.fasta | perl -pe 's/^(.*?[Cc]hr\d+)_/\1\t/; s/_/\t/g' > clust_8.tsv
seqkit seq -n clust9_apipr.MO19963523.gnm1.Chr04_51329060_394_193_seqs.fasta | perl -pe 's/^(.*?[Cc]hr\d+)_/\1\t/; s/_/\t/g' > clust_9.tsv

# Get clusters with the original data as tables to do alignment trimming via ML
for file in clust_[6-9]*.tsv; do # EDIT FOR NUMBER OF CLUSTERS
    ref=`echo $file | perl -pe 's/clust_(.+).tsv/clust_$1_matched.tsv/'`;
    echo $ref
    awk 'BEGIN {FS=OFS="\t"} 
        NR==FNR {key = $1 FS $2 FS $3 FS $4; keys[key] = 1; next} 
        {key = $1 FS $2 FS $3 FS $4} 
        key in keys' $file ../../../ultra.MO19963523.gnm1.814V.p3000.tsv > $ref
done

# After running the ML workflow to output tsv files of the centromeric sequences, we can now make FASTA to look at the alignments and create a consensus sequence
for filename in clust[1-4]_193_all-filters_repeats.tsv; do
    ref=`echo $filename | perl -pe 's/(clust\d+_193_all-filters)_repeats.tsv/$1.fna/'`;
    echo $ref
    awk '{print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $filename > $ref
done

# Delete the first two lines from a batch of files Using sed (edit in-place with backup)
# I could save without a header to optimize this process in the future

for file in clust[0-9]_193_*-filters_*.tsv; do
  sed -i.bak '1,2d' "$file"   # removes lines 1 and 2, keeps a .bak backup
done



### All in One: Get consensus sequences from batch of clustal omega alignments

*Clustal Omega and EMBOSS must be installed*
*This uses the 'setcase' parameter of cons to Set the threshold for the positive matches above which the consensus is is upper-case and below which the consensus is in lower-case.*
*The plurality parameter is set by $threshold. Edit the function according to your needs.*
*You can cut the alignment anywhere using cons -sbegin1 -send1*

    for file in clust[0-9]_193_*-filters_*.fna; do
        echo $file
        seqs=$(grep -c '^>' $file)
        threshold=$(python3 -c "import math; print(max(1, math.ceil($seqs * 0.4)))") # for 'cons' command
        echo $threshold

        ref=`echo $file | perl -pe 's/(.+)\.fna/$1/'`; # EDIT ACCORDINGLY
        #clustalo -i $file -o $ref.aln --outfmt=fa --force && 
        cons -sequence $ref.aln -outseq $ref.cons -plurality $threshold -setcase $threshold # add these parameters if needed, -plurality $threshold -setcase $threshold
    done
    
# revised set of consensus sequences from ML workflow
# Viewing an alignment of these and the previously generated consensus sequences, we can see that the ML workflow refined the consensus sequences with no loss of information.
# In fact, outliers (i.e. sequences that were skewing alignments) can be found in their respective files for each cluster.
>clust1.cons
AAAACTTTTCTTTTTTGGTTATTATTACTTTAATTTTTTTTCCTTACCATTAAATTAACC
TTCTTAACTTTTATTTTACCTTTTTTTAACTTTTTTTTTACTTTGTTTTCTTTTTTTTAT
TTTTTTTTTTAAACCAATTTTACCTTATTAAATTTTTAAAATTTAAGGTTTTGAAAACTT
TTTTTTATTTTTG
>clust2.cons
AAAAAATTTTACCTTTTTTTAACCAAATTATTCCATTTTTCCAAAATCATAAACATTTTT
CCAAACCTTAATTTTTCAAAATTTTATTACGTTAATTTTGTTTTAAAAATTTTTTTAAAA
TTATTAACCAAATTCTTAAATTAGnTTTTATTTTTTCAATTCAAAGTTTAGATTTTTTTT
TTATTTTTTTTTT
>clust3.cons
AAAAAGATTTAAAATATAATGACCGAGAAAGTACAGTTTTACAAAACCATAAAATTGTTT
CCAAAGCCTAAAGTTTCaGAATTTCAGCATACTAAATCTGTATTACAAAATATATTAAAA
ATAGTAAACGAAGTCATAAAATAACTTATTATATGTCnAATTAAAGCTTAGTATGTCTAT
TTTATCTTTATGT
>clust4.cons
AAAATAGACATCTTGAGCTTTAATTTGACATATAATATGTTATTTTATGATTTCGTTTAC
TATTTTTATTATATTTCTTAAGACAGATTTACTACACTATAATTCTGAAAATTTCTGTTG
TGTAATCATTTTTATAATTTTGGAAAATTGTACTACCTTGGTTATAATAATTAAAATATT
TTTACCTAATAAT
>Now from apipr...
>clust6_193_all-filters_repeats.cons
AAAAATAAAATAGACATCTTGAGCTTTAATTTGACATATAATAATCTATTTTATGATTTC
GTTTACTAATTTTAATATATTTATCAAAACAGATTTATTGTATTAGAATTCTGAAAATTT
CTGTTGTGCAACTGTTTTTATAAATTTAGAAAACTGTACTACCTCAGTTATAATAATTTA
AATGTTTTTGTCT
>clust7_193_ML-filters_repeats.cons
AAAAACATTTAAAATATTATGCACGAGAAAGTACAGTTTTATAAAAATGTAAAACCGTTT
CAAGAGCTTAAAGTTTCAGAATTTCAGCATACTAxAATCTGTGTTACTAAATATATTAAA
AATAGTAAAAGAAGTCATAAAATAACTTATTATATGTCAAATTAAAGCTTAGTATGTCTA
TTTTATGTTTATGT
>Removed cluster 8 as it was too small and not informative
>clust9_193_Aln-filters_repeats.cons
AAAACAGATTTATTGTATTAGAATTCTGAAAATTTCTGTTGTGCAACTGATTTTATAATT
TTGGAAAACTGTACTACCTTGTTTATAATAATTTAAATATTTTCACCTAATAATGAAATA
GACATCTTAAGCTTTAATTTGACATATTACAAGTTATTTTATGAATTCGTTTACCATTTT
TATTATATTTCTT

# BLAST searches
  for path in blastdb/*.nsq; do 
    file=`basename $path .nsq`
    echo "WORKING on $file"
    blastn -query query/apios_ALL.cons -db blastdb/$file -out blastout/cent_rpt.x.$file.bln \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" &
  done

# Summarize results
  for filepath in blastout/*bln; do 
    echo $filepath
    cut -f1 $filepath | sort | uniq -c;
    echo
  done

# Find the most frequent target sequences
  cat blastout/cent_rpt.x.Apiam.bln | awk '$3>=90 $4==193 {print $13}' | sort | uniq -c | sort -n | tail -1 | awk '{print $2}'
  # Removing the length ($4) constraint may be necessary to obtain results
    >Apiam_freq
    AAAAATATTTAAATTATTATAACCAAGGTAGTACAGTTTTCCAAAATTATAAAACTG-TTCCAAAGCCTAAAGTTTCAGAATTTCAGCATACTAAATATATGTTAAGAAATATATTAAAAATAGTAAACGAAGTCATAAAATAACTTATTATATGTCGAATTAAATCTTAGTATGTCTATTTTATCTTTATGT
  cat blastout/cent_rpt.x.Apipr.bln | awk '$3>=90 && $4==193 {print $13}' | sort | uniq -c | sort -n | tail -1 | awk '{print $2}'
  # Removing the length ($4) constraint may be necessary to obtain results
    >Apipr_freq
    AAAAATATTTAAATTATTATAACCAAGGTAGTACAGTTTTCCAAAATTATAAAACTG-TTCCAAAGCCTAAAGTTTCAGAATTTCAGCATACTAAATATATGTTAAGAAATATATTAAAAATAGTAAACGAAGTCATAAAATAACTTATTATATGTCGAATTAAATCTTAGTATGTCTATTTTATCTTTATGT
# These are identical sequences, so we can use the same one for both Apiam and Apipr

# Add the most frequents matching sequences to the query and rerun
>Apios_freq
AAAAATATTTAAATTATTATAACCAAGGTAGTACAGTTTTCCAAAATTATAAAACTG-TTCCAAAGCCTAAAGTTTCAGAATTTCAGCATACTAAATATATGTTAAGAAATATATTAAAAATAGTAAACGAAGTCATAAAATAACTTATTATATGTCGAATTAAATCTTAGTATGTCTATTTTATCTTTATGT
>clust1.cons
AAAACTTTTCTTTTTTGGTTATTATTACTTTAATTTTTTTTCCTTACCATTAAATTAACC
TTCTTAACTTTTATTTTACCTTTTTTTAACTTTTTTTTTACTTTGTTTTCTTTTTTTTAT
TTTTTTTTTTAAACCAATTTTACCTTATTAAATTTTTAAAATTTAAGGTTTTGAAAACTT
TTTTTTATTTTTG
>clust2.cons
AAAAAATTTTACCTTTTTTTAACCAAATTATTCCATTTTTCCAAAATCATAAACATTTTT
CCAAACCTTAATTTTTCAAAATTTTATTACGTTAATTTTGTTTTAAAAATTTTTTTAAAA
TTATTAACCAAATTCTTAAATTAGnTTTTATTTTTTCAATTCAAAGTTTAGATTTTTTTT
TTATTTTTTTTTT
>clust3.cons
AAAAAGATTTAAAATATAATGACCGAGAAAGTACAGTTTTACAAAACCATAAAATTGTTT
CCAAAGCCTAAAGTTTCaGAATTTCAGCATACTAAATCTGTATTACAAAATATATTAAAA
ATAGTAAACGAAGTCATAAAATAACTTATTATATGTCnAATTAAAGCTTAGTATGTCTAT
TTTATCTTTATGT
>clust4.cons (may cut this one out)
AAAATAGACATCTTGAGCTTTAATTTGACATATAATATGTTATTTTATGATTTCGTTTAC
TATTTTTATTATATTTCTTAAGACAGATTTACTACACTATAATTCTGAAAATTTCTGTTG
TGTAATCATTTTTATAATTTTGGAAAATTGTACTACCTTGGTTATAATAATTAAAATATT
TTTACCTAATAAT
>clust6_193_all-filters_repeats.cons (may cut this one out)
AAAAATAAAATAGACATCTTGAGCTTTAATTTGACATATAATAATCTATTTTATGATTTC
GTTTACTAATTTTAATATATTTATCAAAACAGATTTATTGTATTAGAATTCTGAAAATTT
CTGTTGTGCAACTGTTTTTATAAATTTAGAAAACTGTACTACCTCAGTTATAATAATTTA
AATGTTTTTGTCT
>clust7_193_ML-filters_repeats.cons (may cut this one out)
AAAAACATTTAAAATATTATGCACGAGAAAGTACAGTTTTATAAAAATGTAAAACCGTTT
CAAGAGCTTAAAGTTTCAGAATTTCAGCATACTAxAATCTGTGTTACTAAATATATTAAA
AATAGTAAAAGAAGTCATAAAATAACTTATTATATGTCAAATTAAAGCTTAGTATGTCTA
TTTTATGTTTATGT
>clust9_193_Aln-filters_repeats.cons
AAAACAGATTTATTGTATTAGAATTCTGAAAATTTCTGTTGTGCAACTGATTTTATAATT
TTGGAAAACTGTACTACCTTGTTTATAATAATTTAAATATTTTCACCTAATAATGAAATA
GACATCTTAAGCTTTAATTTGACATATTACAAGTTATTTTATGAATTCGTTTACCATTTT
TATTATATTTCTT

# After some furthee analysis, I decided to use the following consensus sequences for the final BLAST searches.
# Several consensus sequences were redundant when compared to reverse-complemented sequences, so I removed them.
# The remaining sequences were the most informative and unique.

>Apios_freq
AAAAATATTTAAATTATTATAACCAAGGTAGTACAGTTTTCCAAAATTATAAAACTG-TTCCAAAGCCTAAAGTTTCAGAATTTCAGCATACTAAATATATGTTAAGAAATATATTAAAAATAGTAAACGAAGTCATAAAATAACTTATTATATGTCGAATTAAATCTTAGTATGTCTATTTTATCTTTATGT
>clust6_193
AAAAATAAAATAGACATCTTGAGCTTTAATTTGACATATAATAATCTATTTTATGATTTC
GTTTACTAATTTTAATATATTTATCAAAACAGATTTATTGTATTAGAATTCTGAAAATTT
CTGTTGTGCAACTGTTTTTATAAATTTAGAAAACTGTACTACCTCAGTTATAATAATTTA
AATGTTTTTGTCT
>clust7_193
AAAAACATTTAAAATATTATGCACGAGAAAGTACAGTTTTATAAAAATGTAAAACCGTTT
CAAGAGCTTAAAGTTTCAGAATTTCAGCATACTAxAATCTGTGTTACTAAATATATTAAA
AATAGTAAAAGAAGTCATAAAATAACTTATTATATGTCAAATTAAAGCTTAGTATGTCTA
TTTTATGTTTATGT
>clust9_193
AAAACAGATTTATTGTATTAGAATTCTGAAAATTTCTGTTGTGCAACTGATTTTATAATT
TTGGAAAACTGTACTACCTTGTTTATAATAATTTAAATATTTTCACCTAATAATGAAATA
GACATCTTAAGCTTTAATTTGACATATTACAAGTTATTTTATGAATTCGTTTACCATTTT
TATTATATTTCTT

# So I consulted with the team and we decided to use the following consensus sequences for the final BLAST searches,
# since Apipr_clust6_193, Apipr_clust7_194, Apipr_clust9_193 are very similar — just frame-shifted and/or with some inversion, 
# and Apipr_clust6_193 and Apipr_clust9_193 are 96% identical.

>Apios_freq
AAAAATATTTAAATTATTATAACCAAGGTAGTACAGTTTTCCAAAATTATAAAACTG-TTCCAAAGCCTAAAGTTTCAGAATTTCAGCATACTAAATATATGTTAAGAAATATATTAAAAATAGTAAACGAAGTCATAAAATAACTTATTATATGTCGAATTAAATCTTAGTATGTCTATTTTATCTTTATGT
>clust6_193
AAAAATAAAATAGACATCTTGAGCTTTAATTTGACATATAATAATCTATTTTATGATTTC
GTTTACTAATTTTAATATATTTATCAAAACAGATTTATTGTATTAGAATTCTGAAAATTT
CTGTTGTGCAACTGTTTTTATAAATTTAGAAAACTGTACTACCTCAGTTATAATAATTTA
AATGTTTTTGTCT

# BLAST searches
  for path in blastdb/*.nsq; do 
    file=`basename $path .nsq`
    echo "WORKING on $file"
    blastn -query query/apios_ALL.cons -db blastdb/$file -out blastout/cent_rpt.x.$file.bln \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" &
  done

# Summarize results
  for filepath in blastout/*bln; do 
    echo $filepath
    cut -f1 $filepath | sort | uniq -c;
    echo
  done

blastout/cent_rpt.x.Apaim.bln
   1042 Apiam_freq
    679 Apipr_clust6_193

blastout/cent_rpt.x.Apipr.bln
   2848 Apiam_freq
   3928 Apipr_clust6_193

  
# Extract coords into bed files
  cat blastout/cent_rpt.x.Apiam.bln | 
    awk -v OFS="\t" '$10>$9 {print $2, $9-1, $10}
                     $10<$9 {print $2, $10, $9-1}' | 
    sort -k1,1 -k2n,2n > blastout/cent_rpt.x.Apiam.bed

  cat blastout/cent_rpt.x.Apipr.bln |
    awk -v OFS="\t" '$10>$9 {print $2, $9-1, $10}
                     $10<$9 {print $2, $10, $9-1}' | 
    sort -k1,1 -k2n,2n > blastout/cent_rpt.x.Apipr.bed

# Merge hits within specified distance of one another
  for width in 10 25 50 100; do
    echo "Merge width: ${width}000"
    bedtools merge -i blastout/cent_rpt.x.Apiam.bed -d ${width}000 > blastout/Apiam_merged_${width}k.bed
    bedtools merge -i blastout/cent_rpt.x.Apipr.bed -d ${width}000 > blastout/Apipr_merged_${width}k.bed
  done

  for filepath in blastout/*merged*; do 
    base=`basename $filepath`
    echo $base
    awk -v OFS="\t" 'NR==1 {span=$3-$2; prev=$1; grand+=$3-$2} 
                     NR>1 && $1==prev {span+=$3-$2; grand+=$3-$2} 
                     NR>1 && $1!=prev {print prev, span; prev=$1; span=$3-$2; grand+=$3-$2} 
                     END{print $1, span; print "Total:", grand}' $filepath
    echo
  done

Apiam_merged_100k.bed
apiam.LA2127.gnm1.Chr01 1042988
apiam.LA2127.gnm1.Chr02 1561246
apiam.LA2127.gnm1.Chr03 187413
apiam.LA2127.gnm1.Chr04 1843175
apiam.LA2127.gnm1.Chr05 647006
apiam.LA2127.gnm1.Chr06 275191
apiam.LA2127.gnm1.Chr07 1769626
apiam.LA2127.gnm1.Chr08 601723
apiam.LA2127.gnm1.Chr09 4084525
apiam.LA2127.gnm1.Chr10 1353037
apiam.LA2127.gnm1.Chr11 24756
Total:  13390686

Apiam_merged_10k.bed
apiam.LA2127.gnm1.Chr01 214980
apiam.LA2127.gnm1.Chr02 263892
apiam.LA2127.gnm1.Chr03 73720
apiam.LA2127.gnm1.Chr04 242651
apiam.LA2127.gnm1.Chr05 100015
apiam.LA2127.gnm1.Chr06 96440
apiam.LA2127.gnm1.Chr07 193319
apiam.LA2127.gnm1.Chr08 102417
apiam.LA2127.gnm1.Chr09 339836
apiam.LA2127.gnm1.Chr10 116543
apiam.LA2127.gnm1.Chr11 24756
Total:  1768569

Apiam_merged_25k.bed
apiam.LA2127.gnm1.Chr01 263784
apiam.LA2127.gnm1.Chr02 408287
apiam.LA2127.gnm1.Chr03 93856
apiam.LA2127.gnm1.Chr04 570242
apiam.LA2127.gnm1.Chr05 198704
apiam.LA2127.gnm1.Chr06 142419
apiam.LA2127.gnm1.Chr07 708140
apiam.LA2127.gnm1.Chr08 170379
apiam.LA2127.gnm1.Chr09 651502
apiam.LA2127.gnm1.Chr10 322401
apiam.LA2127.gnm1.Chr11 24756
Total:  3554470

Apiam_merged_50k.bed
apiam.LA2127.gnm1.Chr01 723573
apiam.LA2127.gnm1.Chr02 874454
apiam.LA2127.gnm1.Chr03 122951
apiam.LA2127.gnm1.Chr04 1188850
apiam.LA2127.gnm1.Chr05 236612
apiam.LA2127.gnm1.Chr06 142419
apiam.LA2127.gnm1.Chr07 1166667
apiam.LA2127.gnm1.Chr08 309709
apiam.LA2127.gnm1.Chr09 1254247
apiam.LA2127.gnm1.Chr10 696881
apiam.LA2127.gnm1.Chr11 24756
Total:  6741119

Apipr_merged_100k.bed
apipr.MO19963523.gnm1.Chr01     3213820
apipr.MO19963523.gnm1.Chr02     5069353
apipr.MO19963523.gnm1.Chr03     2920286
apipr.MO19963523.gnm1.Chr04     3090649
apipr.MO19963523.gnm1.Chr05     3041361
apipr.MO19963523.gnm1.Chr06     3727382
apipr.MO19963523.gnm1.Chr07     3998101
apipr.MO19963523.gnm1.Chr08     3394919
apipr.MO19963523.gnm1.Chr09     3079726
apipr.MO19963523.gnm1.Chr10     3491479
apipr.MO19963523.gnm1.Chr11     3490388
Total:  38517464

Apipr_merged_10k.bed
apipr.MO19963523.gnm1.Chr01     890925
apipr.MO19963523.gnm1.Chr02     1038529
apipr.MO19963523.gnm1.Chr03     526939
apipr.MO19963523.gnm1.Chr04     744251
apipr.MO19963523.gnm1.Chr05     634227
apipr.MO19963523.gnm1.Chr06     839518
apipr.MO19963523.gnm1.Chr07     1098201
apipr.MO19963523.gnm1.Chr08     671578
apipr.MO19963523.gnm1.Chr09     575182
apipr.MO19963523.gnm1.Chr10     750838
apipr.MO19963523.gnm1.Chr11     671071
Total:  8441259

Apipr_merged_25k.bed
apipr.MO19963523.gnm1.Chr01     1330238
apipr.MO19963523.gnm1.Chr02     2139771
apipr.MO19963523.gnm1.Chr03     930538
apipr.MO19963523.gnm1.Chr04     1158360
apipr.MO19963523.gnm1.Chr05     948783
apipr.MO19963523.gnm1.Chr06     1383269
apipr.MO19963523.gnm1.Chr07     2071759
apipr.MO19963523.gnm1.Chr08     1143734
apipr.MO19963523.gnm1.Chr09     983947
apipr.MO19963523.gnm1.Chr10     1183197
apipr.MO19963523.gnm1.Chr11     1196150
Total:  14469746

Apipr_merged_50k.bed
apipr.MO19963523.gnm1.Chr01     2250718
apipr.MO19963523.gnm1.Chr02     3500094
apipr.MO19963523.gnm1.Chr03     1586609
apipr.MO19963523.gnm1.Chr04     1550877
apipr.MO19963523.gnm1.Chr05     1619357
apipr.MO19963523.gnm1.Chr06     2046675
apipr.MO19963523.gnm1.Chr07     3006470
apipr.MO19963523.gnm1.Chr08     2180731
apipr.MO19963523.gnm1.Chr09     1776014
apipr.MO19963523.gnm1.Chr10     1956822
apipr.MO19963523.gnm1.Chr11     2231942
Total:  23706309


# Get sequence (e.g. chromosome) sizes from FASTA files
  for file in *.fna; do
    base=`basename $file .fna`
    echo $base
    seqkit fx2tab -n -l $file > ${base}_len.txt
  done




--------------------------
# making jitter plots in R from BLAST output (chromosome lengths files is needed!)
---
setwd("C:/Users/bdjor/Desktop/temp/apios/blastout") # Set your working directory

blast <- read.table("cent_rpt.x.Apiam.bln", header=FALSE, sep="\t",
  col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore","sseq"))

chr_lengths <- read.table("apiam.LA2127.gnm1.KVMK.genome_main_len.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chr_lengths) <- c("sseqid", "length")


# Change chromosome naming 
blast$sseqid <- sub(".*\\.", "", blast$sseqid)
chr_lengths$sseqid <- sub(".*\\.", "", chr_lengths$sseqid)

# merge
blast <- merge(blast, chr_lengths, by = "sseqid")

# Order chromosomes
blast$sseqid <- factor(blast$sseqid, levels = paste0("Chr", sprintf("%02d", 1:11)))


library(ggplot2)
library(dplyr)

# # Get list of chromosomes in your data
# chromosomes <- unique(blast$sseqid)

# # Loop over each chromosome
# for (chr in chromosomes) {
#   chr_data <- filter(blast, sseqid == chr)
  
#   p <- ggplot(chr_data, aes(x = sstart, y = qseqid, color = qseqid)) +
#     geom_jitter(height = 0.3, size = 1.2, alpha = 0.7) +
#     scale_color_viridis_d(option = "C") +  # Use a colorblind-friendly palette
#     labs(
#       title = paste("BLAST Query Start Positions –", chr),
#       x = "Genomic Position",
#       y = "Query ID",
#       color = "Query"
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       axis.text.y = element_text(size = 8),
#       legend.position = "right"
#     )
  
#   # Save each plot as a PDF
#   ggsave(filename = paste0("blast_jitter_", chr, ".pdf"), plot = p,
#          width = 10, height = 6)
# }

# Chatgpt gave me some code to use chromosome lengths to set the x-axis limits, and i fixed it. See below for the final code.

chromosomes <- unique(blast$sseqid)

for (chr in chromosomes) {
  chr_data <- blast[blast$sseqid == chr, ]
  #chr_data <- filter(blast, sseqid == chr)
  chr_length <- unique(chr_data$length.y)  # should be a single number
  
  p <- ggplot(chr_data, aes(x = sstart, y = qseqid, color = qseqid)) +
    geom_jitter(height = 0.3, size = 1.2, alpha = 0.7) +
    scale_color_viridis_d(option = "C") +
    coord_cartesian(xlim = c(0, chr_length)) +
    labs(
      title = paste("BLAST Query Start Positions –", chr),
      x = "Position (bp)",
      y = "Query ID",
      color = "Query"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  ggsave(filename = paste0("blast_jitter_", chr, ".pdf"), plot = p,
         width = 10, height = 6)
}

# Now using plotly to make interactive plots with some added dimension

library(ggplot2)
library(plotly)
library(dplyr)
library(htmlwidgets)

# 
chr_lengths <- read.table("chrom_lengths.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chr_lengths) <- c("sseqid", "chr_length")

# Change chromosome naming 
blast$sseqid <- sub(".*\\.", "", blast$sseqid)
chr_lengths$sseqid <- sub(".*\\.", "", chr_lengths$sseqid)

# merge
blast <- merge(blast, chr_lengths, by = "sseqid")

# Order chromosomes
blast$sseqid <- factor(blast$sseqid, levels = paste0("Chr", sprintf("%02d", 1:11)))

# plotting
chromosomes <- unique(blast$sseqid)

for (chr in chromosomes) {
  chr_data <- blast[blast$sseqid == chr, ]
  #chr_data <- filter(blast, sseqid == chr)
  chr_length <- unique(chr_data$length.y)  # should be a single number
  
  p <- ggplot(chr_data, aes(x = sstart, y = qseqid, color = qseqid)) +
    geom_jitter(height = 0.3, size = 1.5, alpha = 0.7) +
    scale_color_viridis_d(option = "C") +
    coord_cartesian(xlim = c(0, chr_length)) +
    labs(
      title = paste("BLAST Query Start Positions –", chr),
      x = "Position (bp)",
      y = "Query ID",
      color = "Query"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )
  
  interactive_plot <- ggplotly(p, tooltip = c("x", "y", "color"))

  # Save as standalone HTML
  saveWidget(interactive_plot, file = paste0("plotly_", chr, ".html"), selfcontained = TRUE)
}

#---- facet wrapped multiplot of each chromosome ----
# FIRST IN THE SHELL, run the following to source your .env file and create the env.r file for R to read
vim env.sh # edit this file to set your environment variables
source env.sh
cat env.sh | sed 's/^export //' > .env_r

# IN R: Set working directory
#setwd("C:/Users/bdjor/Desktop/temp/apios/blastout") # Set your working directory
setwd(apios_path)

# Load variables from .env file
readRenviron(".env_r")
apios_path <- Sys.getenv("PWD")
apiam_bln_path <- Sys.getenv("APIAM_BLAST")
apipr_bln_path <- Sys.getenv("APIPR_BLAST")
apiam_len <- Sys.getenv("APIAM_LEN")
apipr_len <- Sys.getenv("APIPR_LEN")



# Read BLAST output
blast <- read.table("cent_rpt.x.Apiam.bln", header=FALSE, sep="\t",
  col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore","sseq"))

# Read chromosome lengths   
chr_lengths <- read.table("apiam.LA2127.gnm1.KVMK.genome_main_len.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chr_lengths) <- c("sseqid", "chr_length")

# Change chromosome naming 
blast$sseqid <- sub(".*\\.", "", blast$sseqid)
chr_lengths$sseqid <- sub(".*\\.", "", chr_lengths$sseqid)


# Create named vector of colors
# Get unique query IDs
queries <- levels(factor(blast$qseqid))
# Create a color palette
isu_colors <- c("#C8102E", "#F1BE48", "#333333", "#A7A9AC", "#FDB827", "#8C1D40")
# Repeat or truncate to fit number of queries
color_palette <- rep(isu_colors, length.out = length(queries))
names(color_palette) <- queries


# merge
blast <- merge(blast, chr_lengths, by = "sseqid")

max_chr_length <- max(blast$chr_length, na.rm = TRUE)

# Order chromosomes
blast$sseqid <- factor(blast$sseqid, levels = paste0("Chr", sprintf("%02d", 1:11)))

# Prompt for alias for naming and labeling
alias <- readline("Enter a short alias for labeling and output filename: ")

library(ggplot2)
library(patchwork)

p <- ggplot(blast, aes(x = sstart, y = qseqid, color = qseqid)) +
  geom_jitter(height = 0.3, size = 0.8, alpha = 0.7) +

  # This makes each facet respect the chromosome length
  geom_blank(aes(x = chr_length)) +

  #scale_color_viridis_d(option = "C") + # original color options
  scale_color_manual(values = color_palette) +  #  <== apply ISU colors here
  #facet_wrap(~ sseqid, scales = "free_x", ncol = 2, drop = FALSE) + # with variable chr_lengths
  facet_wrap(~ sseqid, scales = "fixed", ncol = 1, drop = FALSE) + # using max length for all facets
  coord_cartesian(xlim = c(0, max_chr_length)) + # fix axis limits

  labs(
    title = paste(alias," BLAST Query Positions"),
    x = "Position (bp)",
    y = "Query ID",
    color = "Query"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 10),
    legend.position = "none"
  )

# Construct output filename
plotout <- paste0("test_blast_faceted_by_chr_length_", alias, ".pdf")
# Save the plot
ggsave(plotout, plot = p, width = 12, height = 16)


# hack to clear the RStudio console, as R has no clear console command

cat("\014")  # This sends a form feed character to the console, which clears it in RStudio

# Side-by-side vertical comparison of Apiam and Apipr BLAST results
# Read BLAST output
blast1 <- read.table(apiam_bln_path, header=FALSE, sep="\t",
  col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore","sseq"))
blast2 <- read.table(apipr_bln_path, header=FALSE, sep="\t",
  col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore","sseq"))

# Read chromosome lengths   
chr_lengths_1 <- read.table(apiam_len, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chr_lengths_1) <- c("sseqid", "chr_length")
chr_lengths_2 <- read.table(apipr_len, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chr_lengths_2) <- c("sseqid", "chr_length")

# Change chromosome naming 
blast1$sseqid <- sub(".*\\.", "", blast1$sseqid)
chr_lengths_1$sseqid <- sub(".*\\.", "", chr_lengths_1$sseqid)
blast2$sseqid <- sub(".*\\.", "", blast2$sseqid)
chr_lengths_2$sseqid <- sub(".*\\.", "", chr_lengths_2$sseqid)


# merge
blast1 <- merge(blast1, chr_lengths_1, by = "sseqid")
blast2 <- merge(blast2, chr_lengths_2, by = "sseqid")

# Get maximum chromosome length for setting x-axis limits
max_chr_length <- max(c(blast1$chr_length, blast2$chr_length), na.rm = TRUE)

# Order chromosomes
chrom_order <- paste0("Chr", sprintf("%02d", 1:11))
blast1$sseqid <- factor(blast1$sseqid, levels = chrom_order)
blast2$sseqid <- factor(blast2$sseqid, levels = chrom_order)

# Load necessary libraries
library(ggplot2)
library(patchwork)

# Define a helper plotting function
plot_blast <- function(df, fill_color = "#C8102E", max_chr_length, individual_title = NULL, show_title = FALSE) {
  ggplot(df, aes(x = sstart)) +
    geom_jitter(aes(y = 0), height = 0.3, size = 0.8, alpha = 0.7, color = fill_color) +
    geom_blank(aes(x = chr_length)) +
    facet_wrap(~ sseqid, scales = "fixed", ncol = 1, drop = FALSE) +
    coord_cartesian(xlim = c(0, max_chr_length)) +
    labs(
      x = NULL,
      y = "193 bp",
      title = if (show_title) individual_title else NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),    # no tick labels
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(angle = 90, size = 10),
      strip.text = element_text(size = 10),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
}

# Define alternating ISU colors (Cardinal and Gold)
isu_colors <- c("#C8102E", "#F1BE48")

# Use color 1 for blast1, color 2 for blast2
p1 <- plot_blast(blast1, fill_color = isu_colors[1], max_chr_length, individual_title = "Apiam", show_title = TRUE)
p2 <- plot_blast(blast2, fill_color = isu_colors[2], max_chr_length, individual_title = "Apios", show_title = TRUE)

# Combine them
final_plot <- p1 | p2

# OR add global title
final_plot <- (p1 | p2) +
  patchwork::plot_annotation(
    title = "Locations of Putative 193 bp Centromeric Repeat",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

# Save
ggsave("patchwork_dual_blasts_ISU.pdf", plot = final_plot, width = 14, height = 20)
