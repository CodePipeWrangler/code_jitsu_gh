# Generate a consensus sequence to locate centromeric regions

Many questions came to mind after finding putative centromeric regions and their respective DNA repeats, such as are the repeats alike? How much do they vary? from location to location in the genome? on the same chromosome? Hence, I pooled putative centromeric repeats and did phylogenetic analyses to get some answers. 

Centromeric repeats in the same species can somewhat vary in sequence as you progress away from the centromere and at neo-centromeres. Furthermore, sometimes you may be looking at the reverse complement of a particular DNA repeat.C hromosomal restructuring throughout evolution can invert DNA sequences or you could be looking at the other strand. Therefore, pooling all repeats of even a single centromeric repeat class (e.g. 92 bp repeats) can skew the alignment significantly and prevent generation of a high quality consensus.

## To illustrate,

### Extract data for specific repeat sizes and create FASTA files

file = ULTRA_FILE_CONVERTED_2_TSV # see [characterize_glycine_centromeres](https://github.com/CodePipeWrangler/code_jitsu_gh/blob/main/genomics/centromeric_characterization/characterize_glycine_centromeres.md) for instructions on how to generate this file

```shell
    awk '$4==91 {print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $file >> my_cent_seqs_91.fn 
    awk '$4==92 {print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $file >> my_cent_seqs_92.fn 
```

*Sample of Alignment of all putative and neo-centromeric 91 bp monomers from Wm82.gnm6 genome assembly.*

*Created via [Clustal Omega](http://www.clustal.org/omega/) (available at the [EMBL_EBI Job Dispatcher](https://europepmc.org/article/MED/38597606) ) *
![glyma gnm6 91 ALL aln](https://github.com/user-attachments/assets/0cfd3318-beab-4061-9ce8-ce2e25be7aa8)*


Short alignments of tandem repeat monomers can be extremely misleading, given ULTRA generates edge artifacts caused by it's arbitrary cut points when it defines what is a repeat monomer. To assess whether real differencs are being seen in the alignment, I simulated a longer tandem array of every repeat (see code below) and realligned them.

     filename=my_cent_seqs_92.fn
     awk '/^>/ {h=$0; print h} !/^>/ {for(i=0; i<4; i++) print $0}' "$filename" > "${filename%.fn}.4X.fn"

Nevertheless, further analysis benefits from extracting repeats from bins of interest and analyze those. For instance, I may pick 1 Mbp around the peak of 91 bp repeats and pool sequences from every chromosome to determine if they are conserved in the genome. Although centromeric sequences can be highly variable between species, it is expected that variation is less than 5-10% within a species. This is evident in *Glycine max*, where the 91 and 92 bp repeats from Wm82, Lee, and FiskebyIII are nearly identical.

### Extract repeat class from ULTRA file and create FASTA per chromosome.

```shell
    x='91' 
    for filename in ULTRA*tsv; do
        ref=`echo $filename | perl -pe 's/ultra\.(.+)\.p3000\.tsv/$1/'` # REFERENCE TO NAME OUTPUT; EDIT PERL REGEX ACCORDING TO YOUR NEEDS
        echo -e "\n$ref";
        for i in {01..20}; do
            maxBin=`awk -v CUT="$x" ''/'[Cc]hr'$i'/ && $4==CUT {print int($2/1000000)}' $filename | sort -n | uniq -c | awk -v max=0 '{if($1>max){want=$2;     max=$1}}END{print want}'`;
            pos1=$(($maxBin*1000000)); pos2=$((pos1+1000000));
            echo -e "Chr$i\t$x\t$pos1\t$pos2";
            awk -v CUT="$x" -v START="$pos1" -v STOP="$pos2" '$4==CUT && $2>=START && $2<=STOP && '/'[Cc]hr'$i'/ {print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $filename >> $ref.$x.chr$i.cent_$maxBin\M.fn ;
        done
    done
```


To create a consensus sequence to be used for predicting centromeric regions of the genome, I extracted putative centromeric repeats and analyzed the sequences further. Below I detail two complementary approaches to obtain a high-quality consensus:

### Cluster sequences and extract representatives for centromeric repeats

I took these FASTA collections of centromeric repeats and used the program [mmseqs2](https://github.com/soedinglab/MMseqs2) to cluster sequences and extract representative sequences for the two types

  Create mmseqs database from FASTA file for clustering

    mmseqs createdb input.fasta inputDB

  Cluster using mode 1 (It is good to try the other cluster modes 0-3 to compare results)

    mmseqs cluster inputDB clusteredDB tmp cluster-mode 1 --threads 8

  Extract representative sequences

    mmseqs createtsv inputDB inputDB clusteredDB clustered_pairs.tsv
    mmseqs createsubdb clusteredDB inputDB repsDB
    mmseqs convert2fasta repsDB reps.fasta

  See how many clusters were formed
  
  ```shell
      grep '>' reps.fasta | wc -l
  ```
    
  See representative sequence IDs for clusters

  ```shell
    awk '{print $1}' clustered_pairs.tsv | uniq -c | sort -k2,2
  ```

  Extract sequences clustering with a representative fasta

```shell
    rm tmp_ids.txt
    REP="glyma.Wm82.gnm6.Gm08_30848187_11811_92" # INSERT SEQID HERE
    awk -v rep="$REP" '$1==rep {print $2}' clustered_pairs.tsv > tmp_ids.txt
    echo "$REP" >> tmp_ids.txt
    seqkit grep -f tmp_ids.txt input.fasta > clust_${REP}_seqs.fasta # INSERT PATH 2 ORIGINAL FILE CREATING DB HERE
```
## 1. Manual Trimming of Alignments

Alignments of tandemly repeated monomers can be manually trimmed and curated. This approach relies on inspecting the aligned repeat arrays and selecting regions with high similarity to construct the consensus. It is especially useful when representative monomers are highly diverged or structural irregularities are present.
For this approach it is recommended to pool bins of putative centromeric repeats from 1 Mb regions suspected to be centromeric and simulate arrays to facilitate processing. Once an alignment of these sequences is obtained, find a point within where the alignment looks best and cut it to your target length.

Check back soon for a tutorial on manual trimming of alignments for creating centromeric consensus sequences.

## 2. Machine Learning-Based Trimming

In an alternative approach, I developed a [machine learning workflow](https://github.com/CodePipeWrangler/code_jitsu_gh/blob/main/genomics/centromeric_characterization/aln_trim_by_ML_feature_class.py) that uses a representative centromeric sequence to guide trimming of alignments. The method efficiently detects and retains alignment regions that best represent the underlying monomer structure. This enables rapid generation of consensus sequences from large-scale repeat datasets, bypassing the need for extensive manual curation.

This ML-based strategy offers substantial speed advantages and is scalable for use in other species. While it can stand alone, its outputs can also be refined with manual trimming steps, if needed.


## Published sequences for CentG1-1 and -2 vs. my consensus sequences 

After creating consensus sequences from my alignments of centromeric repeats I checked their similarity to those of the published sequences, CentGm-1 and -2. You can see that the consensus and the published sequences are almost identical besides the edge artifacts.

    CentGm-1
    TGTGAAAAGTTATGACCATTTGAATTTCTCGAGAGCTTCCGTTGTTCAATTTCGAGCGTCTCGATATATTATGCGCCTGAATCGGACATCCG

    CentGm-2
    AGTCAAAAGTTATTGTCGTTTGACTTTTCTCAGAGCTTCCGTTTTCAATTACGAGCGTCTCGATATATTACGGGACTCAATCGGACATCCG

    Glyma_91bp_rep
    AAAAAGTTATTGTCGTTTGAATTTGCTCAGAGCTTCAACATTCAATTTCGAGCGTCTCGATATATTACGGGACTCAATCAGACATCCGAGT

    Glyma_92bp_rep
    AAAAGTTATGACCATTCGAATTTCTCGAGAGCTTCCGTTGTTCAATTTCGAGCGTCTCGATATATTATGTCCCCGAATCGGACATTCGTGTG

![Screenshot 2025-06-15 130714](https://github.com/user-attachments/assets/710ed151-7a07-4541-bcf7-64bbb34563b1)
![Screenshot 2025-06-15 123614](https://github.com/user-attachments/assets/7b93831d-a596-49bc-8110-02fc7e5d9e7d)


--------

* INSERT blast results for both repeat types: where do the two types of repeats localize to? same chromosomes? same copy number?*


*two types of repeats in Glycine max.*

CentGm-1 and 2 show large variation in copy number. On most chromosomes CentGm-1 is dominant, while others show CentGm-2 or codominance.



