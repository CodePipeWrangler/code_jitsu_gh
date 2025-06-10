###### Author: Brandon D. Jordan
###### Started 07/08/2024
###### Obj: Bioinformatics tools from the commandline for centromeric repeat research

*I began this project using the software Tandem Repeat Finder (TRF) to identify tandem repeats in the six perennial glycine species. In this analysis, the program ULTRA for tandem repeat identification replaced TRF because of its promoted utility and speed over TRF. 


### *Many of the commands below are written as one-liners. Individual commands reside between semi-colons (;).*
### *It may be useful while editing commands to separate them on individual lines.*

### Set Environmental variables

file=
ref=`echo $file | perl -pe 's/ultra\.(.+)\.p\d.+\.tsv/$1/'` # Create a shorter reference for the filename (modify regex ad hoc)
#ref=`echo $file | perl -pe 's/ultra\.(.+)\.p3000\.tsv/$1/'`
#*NOTE, some genomes annotate chromosomes with identifiers other than Chr, such as chr, Gm, etc. Check the file format before proceeding with analyses

## Project-specific File management
### Transfers files between local/remote locations via scp

    cd WORKING_DIR  
    #scp FILENAME jordan25@perlmutter.nersc.gov:'/pscratch/sd/j/jordan25/genomics/apios/query'
    scp FILENAME REMOTE@SERVER.LOCATION.EDU:'/FILEPATH/FILE'

### Clean up directories and manage data (02-10-2025)
#### edit filename pre- and suffix at once. Here  I remove 'BV-YZ2020.gnm2.17QQ.' and 'align' from the prefix and suffix respectively

    for file in *align; do 
        newfile="${file#BV-YZ2020.gnm2.17QQ.}" # INSERT PREFIX PATTERN HERE AFTER '#' TO REMOVE
        newfile="${newfile%align}aln" # INSERT SUFFIX PATTERN HERE AFTER '%' TO REMOVE
        mv $file $newfile; 
    done

## Data wrangling
### Use command substitution to initialize a loop. Here I use the $ref variable to provide input files to my loop

    for i in `cut -f1 ../chrom_len/"$ref"_chrom_lengths.txt`; do
        CODEBLOCK
    done

### Subset data and print specific cols with delimiter

    awk '$3>=5000 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $9}' $file

## Working with FASTA files
### Remove new lines from sequences in FASTA

    awk '/^>/ {print (NR==1 ? "" : "\n") $0; next} {printf "%s", $0}' $file 

#### Convert multi-line FASTA to single-lined

    awk '/^>/ {printf("%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <  multiline.fasta >> singleline.fasta

### tandemly duplicate FASTA sequences 4X times and save to new file

    filename=input.fna
    awk '/^>/ {h=$0; print h} !/^>/ {for(i=0; i<4; i++) print $0}' "$filename" > "${filename%.fna}.4X.fna"


#### Seqkit utilities
### Remove duplicate entries from a FASTA using seqkit

    seqkit rmdup -s < seqs.fasta > result.fasta

### Search with list of sequences IDs

    seqkit grep -f id.txt seqs.fasta -o result.fasta

### Search with pattern

    seqkit grep -r -p "pattern_to_exclude" -v seqs.fasta > results.fasta

### Exclude sequences by ID

    seqkit grep -f id.txt -v seqs.fasta > results.fasta

### Subset all FASTA sequences by position

    seqkit subseq -r START_POS:END_POS $file > subset.X1-X2.fasta

## DATA ANALYSIS
#### ON ULTRA RESULTS
### Check the chromosomal naming scheme in ULTRA results

    cut -f1 $file | head

### Check the number of chromosomes represented in results, ignoring scaffolds (scaffold sequences are often denoted by 'sc' or 'scaffold'. 
### However it is always a best practice to check for your case)

    cut -f1 $file | uniq -c /
        | grep -v 'sc' # to remove scaffolds from the list

### Loop to sort ultra.tsv files by cols 1>3 for downstream ease-of-use, and also remove extraneous extensions from upstream. This creates a separate file prefixed with 'tmp'

    for filename in *tsv; do 
        echo; ref=`echo $filename | perl -pe 's/ultra\.(.+)\.p3000\.tsv/$1/'`
        cat $filename | sort -t$'\t' -k1,1 -k3,3n > tmp.$ref.tsv
        echo tmp.$ref.tsv;  
    done

### Histogram of the period distribution (adjust the trailing awk and perl commands respectively for filtering and noise cancellation)

    awk '{print $4}' $file | sort -n | uniq -c | awk '$1>=10 && $2>=60 {print $1 "\t" $2}' | perl -lane 'print $F[1], "\t", "." x int($F[0]/100)'

### Histogram by Chr of the period dist. (adjust the trailing awk and perl commands respectively for filtering and noise cancellation)

    for i in {01..10}; do echo Chr$i ; awk ''/'.Chr'$i'/ {print $4}' $file | sort -n | uniq -c | awk '$1>=10 && $2>=60 {print $1 "\t" $2}' | perl -lane 'print $F[1], "\t", "." x int($F[0]/100)'; echo; done

### Histogram by Chr of a particular repeat dist. by location (adjust the trailing awk and perl commands respectively for filtering and noise cancellation)

    for i in {01..10}; do echo Chr$i ; awk ''/'.Chr'$i'/ && $4==104 {print int($2/1000000)}' $file | sort -n | uniq -c | awk '$1>=1 && $2>=1 {print $1 "\t" $2}' | perl -lane 'print $F[1], "\t", "." x int($F[0]/100)'; echo -e '\n'; done

### Histogram by chr. of satellite array sizes for repeats of a specific length or period range

    for i in {01..10}; do echo Chr$i; awk ''/'Chr'$i'/ && $4>=155 && $4<=156 {print int($3/1000)}' $file | sort -n | uniq -c | awk '$1>=1 && $2>=1 {print $1 "\t" $2}' | perl -lane 'print $F[1], "\t", "." x int($F[0]/10)'; echo; done


### This will print the largest bin (size: 1Mbp) for each Chr that contains the repeat of interest

    x=156
    echo $ref
    for i in {1..7}; do echo Chr$i ;
        awk ''/'Chr'$i'/ && $4=='$x' {print int($2/1000000)}' $file | sort -n | uniq -c | awk -v max=0 '{if($1>max){want=$2; max=$1}}END{print want}';
        echo;
    done

### Extract data for 200 of the largest tandem repeat arrays from each Chr (note these tmp.*tsv files have already been sorted by chr>arraySize as discussed above in this section).
    
    for i in {01..07}; do echo Chr$i ; awk ''/'.Chr'$i'/ {print}' tmp.$file | tail -200 >> sub.$ref.tsv ; done

## Extract sequences for select repeats into FASTA
### Extract 10 of the largest tandem repeat arrays by chr of a given period size into FASTA

    x=156; 
    for i in {01..14}; do echo Chr$i ; awk ''/'.Chr'$i'/ && $4=='$x' {print}' tmp.$ref.tsv | tail -10 | awk '{print ">"$1"_"$2"_"$3"_"$4"\n"$9}' >> sub.$ref.fasta ; done

### Extract data for specific localized repeat sizes and create FASTA

    x='156'; 
    pos1=48000000 ; pos2=52000000
    i=1
    for sub in `echo $x`; do
        awk -v CUT="$sub" -v START="$pos1" -v STOP="$pos2" '$4==CUT && $2>=START && $2<=STOP && '/'Chr'$i'/ {print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $filename >> $ref.$x.chr$i.fasta ;
    done

### This creates FASTA of EACH CHR containing the proposed centromeric regions (largest bin of 'x' bp repeats) for 1+ species (i.e. inputs)
### It also prints the region selected from per chr
### Here I use nested loops to target several genomes

    x='91' 
    for filename in ULTRA*tsv; do
        ref=`echo $filename | perl -pe 's/ultra\.(.+)\.p3000\.tsv/$1/'` # EDIT ACCORDING TO YOUR FILENAME
        echo -e "\n$ref";
        for i in {01..20}; do
            maxBin=`awk -v CUT="$x" ''/'[Cc]hr'$i'/ && $4==CUT {print int($2/1000000)}' $filename | sort -n | uniq -c | awk -v max=0 '{if($1>max){want=$2; max=$1}}END{print want}'`;
            pos1=$(($maxBin*1000000)); pos2=$((pos1+1000000));
            echo -e "Chr$i\t$x\t$pos1\t$pos2";
            awk -v CUT="$x" -v START="$pos1" -v STOP="$pos2" '$4==CUT && $2>=START && $2<=STOP && '/'[Cc]hr'$i'/ {print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $filename >> $ref.$x.chr$i.cent_$maxBin\M.fn ;
        done
    done

### Create single FASTA file from batch of single-sequence files without headers

    input_dir="path_to_your_directory" # Directory containing your DNA sequence files
    output_file="combined_sequences.fasta"

    rm "$output_file" # Create or clear the output file, in case it already exists

    for file in *ext; do 
        filename=$(basename "$file")
        identifier="${filename%.*}"
        echo ">$identifier" >> "$output_file"
        cat "$file" >> "$output_file"
    done

## Data analysis
### [DEPRECATED] Generate consensus sequences for each FASTA (try higher thresholds above the default (0.1))
### This was a script I created to get consensus sequences, and it is available in my github repo to try out.

    for file in *fn; do
        ref=`echo $file | perl -pe 's/(.+)\.fn/$1/'`; echo $ref
        python ~/code_jitsu_gh/genomics/consensus_seq/consensus.py $file 0.1 > $ref.cons
    done

### Generate clustalo alignments for batch of FASTA

    for file in *fna ; do ref=`echo $file | perl -pe 's/(.+)\.fna/$1/'`; echo "\n$ref" ; clustalo -i $file -o $ref.clust.align --outfmt=clu; done

### Generate famsa alignments for batch of FASTA

    for file in *fna ; do ref=`echo $file | perl -pe 's/(.+)\.fna/$1/'`; echo "\n$ref" ; famsa -t 2 $file $file.famsa.aln; done

### All in One: Get consensus sequences from batch of clustal omega alignments

*Clustal Omega and EMBOSS must be installed*
*This uses the 'setcase' parameter of cons to Set the threshold for the positive matches above which the consensus is is upper-case and below which the consensus is in lower-case.*
*The plurality parameter is set by $threshold. Edit the function according to your needs.*
*You can cut the alignment anywhere using cons -sbegin1 -send1*

    for file in *fna; do
        seqs=$(grep -c '^>' $file)
        threshold=$(python3 -c "import math; print(max(1, math.ceil($seqs * 0.2)))")
        echo $threshold

        ref=`echo $file | perl -pe 's/(.+)\..+\.fna/$1/'`; echo "$ref"; # EDIT ACCORDINGLY
        clustalo -i $file -o Type_$ref.aln --outfmt=fa --force && 
        cons -sequence $ref.aln -outseq $ref.cons -plurality $threshold -setcase $threshold
    done

### Get pairwaise comparions with clustalo

    clustalo -i input.fasta --distmat-out=distances.txt --full --percent-id --force

#### On mmseqs files
### Perform cluster analysis of a FASTA database/file


------------------------------------------------------------------------------------------------------------------------
UPDATE 06-10-2025 - 5:08 PM -> SCANNED AND STRUCTURED ENTIRE DOCUMENT, AND I WILL START BACK POLISHING FROM HERE
------------------------------------------------------------------------------------------------------------------------

	cd [mmseqs.directory]
    file=input.fasta
    mode=2

	mkdir clust_mod0 clust_mod1 clust_mod2 clust_mod3 clust_hash.90 db aln prefilter
	mmseqs createdb *fna db/queryDB
	mmseqs cluster db/queryDB clust_mod2/DB_clu tmp --cluster-mode 2
	mmseqs createtsv db/queryDB db/queryDB clust_mod2/DB_clu clust_mod2/DB_clu.tsv
	mmseqs createseqfiledb db/queryDB clust_mod2/DB_clu clust_mod2/DB_clu_seq
	mmseqs result2flat db/queryDB db/queryDB clust_mod2/DB_clu_seq clust_mod2/DB_clu_seq.fasta
	mmseqs createsubdb clust_mod2/DB_clu db/queryDB clust_mod2/DB_clu_rep         
	mmseqs convert2fasta clust_mod2/DB_clu_rep clust_mod2/DB_clu_rep.fasta

#### one-liner to run all clustering modes ( after setting up directories) at once except clust_hash
	
	mkdir clust_mod0 clust_mod1 clust_mod2 clust_mod3 clust_hash.90 db aln prefilter;
	for i in {0..3}; do # EDIT THE CLUSTERING MODES OUTPUT HERE
	mmseqs createdb *fnn db/queryDB && mmseqs cluster db/queryDB clust_mod$i/DB_clu tmp --cluster-mode $i && mmseqs createtsv db/queryDB db/queryDB clust_mod$i/DB_clu clust_mod$i/DB_clu.tsv && mmseqs createseqfiledb db/queryDB clust_mod$i/DB_clu clust_mod$i/DB_clu_seq && mmseqs result2flat db/queryDB db/queryDB clust_mod$i/DB_clu_seq clust_mod$i/DB_clu_seq.fasta && mmseqs createsubdb clust_mod$i/DB_clu db/queryDB clust_mod$i/DB_clu_rep  && mmseqs convert2fasta clust_mod$i/DB_clu_rep clust_mod$i/DB_clu_rep.fasta ;done

### check number of clusters in each representative sequence file for clusterings

    for f in clust_*; do echo $f; grep '>' $f/DB_clu_rep.fasta | wc -l; done

### find which clusters are the largest (all representartive sequences of clusters in DB_clu_rep.fasta can be used to search for homology agaisnt the genome of interest)

	awk '{print $1}' clust*0/*tsv | uniq -c | sort -r | head -30

### find lines in the cluster file where the end of the ID contains '_7'

	awk '$1~/.+_[7]$/ {print $1"\t"$2}' *tsv



### This creates a table of the proposed centromeric regions (largest bin of 156 bp repeats) for 1+ species (i.e. inputs)

    for file in *YOURFILE*tsv; do         
    ref=`echo $file | perl -pe 's/ultra\.Zm-(.+)-REFERENCE-NAM-1.0\.p3000\.tsv/$1/'`;echo "loc\tlen\tbin_start\tbin_end";
    echo "\nloc\tlen\tbin_start\tbin_end";
    for i in {1..10}; do
    maxBin=`awk ''/'chr'$i'/ && $4==156 {print int($2/1000000)}' $file | sort -n | uniq -c | awk -v max=0 '{if($1>max){want=$2; max=$1}}END{print want}'`;
    x='156'; pos1=$(($maxBin*1000000)); pos2=$((pos1+1000000));
    echo "Chr$i\t$x\t$pos1\t$pos2";
    done;
    done

## On BLAST results
### Run a BLASTn search (on Ceres)

    ml blast+ # load the module
    makeblastdb -in cerca.ISC453364.gnm3.GWXB.genome_main.fna -parse_seqids -dbtype nucl -out blastdb/cerca.ISC453364.gnm3.GWXB.genome_main.fna

    blastn -db blastdb/ -query cons*866*.fnn -task blastn -outfmt "7" > blast.cercavs866.txt

### Analyze blast data faster

    i=01; echo Chr$i ; awk '$2~'/'.Chr'$i'/ && $3>=80 && $4>=100 && $9>=4000000 && $10<=9000000 {print int($9/1000000)}' $file | sort -n | wc -l

    for i in {01..07}; do echo Chr$i ; awk '$2~'/'.Chr'$i'/ && $3>=75 && $4>=100 {print int($9/1000000)}' $file | sort -n | uniq -c | awk '$1>=1 && $2>=1 {print $1 "\t" $2}' | perl -lane 'print $F[1], "\t", "." x int($F[0]/1000)'; done
    
    file=blast.cerchi-Cerciscons.txt 

#### Get counts within a region

    - i=01; l=11000000; r=17000000;echo Chr$i ; awk -v LEFT="$l" -v RIGHT="$r" '$2~'/'.Chr'$i'/ && $3>=80 && $4>=100 && $9>=LEFT && $10<=RIGHT {print int($9/1000000)}' $file | sort -n | wc -l ; awk -v LEFT="$l" -v RIGHT="$r" '$2~'/'.Chr'$i'/ && $3>=80 && $4>=100 && $9>=LEFT && $10<=RIGHT {print}' $file > test.Chr$i.txt
