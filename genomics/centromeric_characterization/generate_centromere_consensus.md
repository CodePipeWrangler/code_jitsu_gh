# Generate a consensus sequence to locate centromeric regions

To create a consensus sequence to be used for predicting centromeric regions of the genome, I extracted putative centromeric repeats and analyzed the sequences further.

Are the repeats alike? transposons? How much do they vary? from location to location in the genome? on the same chromosome? Many questions can come to mind, so I grabbed all of the results for putative centromeric repeats and did phylogenetic analyses to get some answers. For *Glycine max* these repeats are similiar enough from place to place that we can pool them all together and still proceed with the analysis without issue. However, there are many cases where it benefits clarity in the results to take repeats from bins of interest and then analyze those. For instance, I may pick 1 Mbp around the peak of 91 bp repeats and pool sequences from every chromosome to determine if they are conserved in the genome. Although centromeric sequences can be highly variable between species, it is expected that variation is less than 5-10% within a species. This is evident in *Glycine max*, where the 91 and 92 bp repeats from Wm82, Lee, and FiskebyIII are nearly identical.

### Extract data for specific repeat sizes and create FASTA files
    awk '$4==91 {print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $file >> my_cent_seqs_91.fn 
    awk '$4==92 {print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $file >> my_cent_seqs_92.fn 

*cluster results and alignment of many 91-92 bp sequences.*
*alignments of putative and neo-centromeric 91-92 Bp repeat monomers.*
    - 91 alignments 
![glyma gnm6 91 ALL aln](https://github.com/user-attachments/assets/638439fd-60dd-4550-863a-a4346b6ea0b8)

Short alignments of tandem repeat monomers can be misleading, given ULTRA defines what is a repeat monomer somewhat arbitrarily. To assess whether real differencs are beign seen in the alignment, I increased the size of every repeat by four times (see code below) and realligned them.

     filename=my_cent_seqs_92.fn
     awk '/^>/ {h=$0; print h} !/^>/ {for(i=0; i<4; i++) print $0}' "$filename" > "${filename%.fn}.4X.fn"

<img width="804" alt="glyma gnm6 91 ALL 4X aln" src="https://github.com/user-attachments/assets/d14fbb7d-1f2a-499b-b297-587a3e3a5428" />

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

    for x in 91 92; do
      grep '>' results/$x.wgd.4X_reps.fasta | wc -l
    done

  See representative sequence IDs for clusters

    for x in 91 92; do awk '{print $1}' clusters/$x*tsv | uniq -c | sort -k2,2 ; done

  Extract sequences clustering with a representative fasta

    rm tmp_ids.txt
    REP="glyma.Wm82.gnm6.Gm08_30848187_11811_92" # INSERT SEQID HERE
    awk -v rep="$REP" '$1==rep {print $2}' clusters/92_clustered_pairs.tsv > tmp_ids.txt
    echo "$REP" >> tmp_ids.txt
    seqkit grep -f tmp_ids.txt glyma.Wm82.gnm6.JFPQ.92.wgd.4X.fn > clust_${REP}_seqs.fasta # INSERT PATH 2 ORIGINAL FILE CREATING DB HERE

### Published sequences for CentG1-1 and -2 vs. my consensus sequences 

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

You can see that the representative sequences and the published CentGm-1 and -2 sequences are almost identical besides the shift in their cut point for the monomer. This is due to how a tandem repeat is cut to define a monomer in each respective study. In this study, ULTRA defined repeat periods.

Sometimes you may see a repeat type on a chromosome or at a neo-centromere that looks to be significantly different in composition from others in the set. It could be there are more than one repeat type in the data, even possibly of the same size. However, chromosomal restructuring throughout evolution can invert DNA sequences or you could be looking at the other strand. To investigate this we can take the reverse complements of the consensus sequences and check their similarity. In this case, the reverse complements are 91% and 92% similar to the 91 bp and 92 bp repeats respectively.

---

* INSERT blast results for both repeat types: where do the two types of repeats localize to? same chromosomes? same copy number?*


*two types of repeats in Glycine max.*

CentGm-1 and 2 show large variation in copy number. On most chromosomes CentGm-1 is dominant, while others show CentGm-2 or codominance.



