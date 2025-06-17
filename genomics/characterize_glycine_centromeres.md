# Characterization of Centromeric Repeats in Glycine Max and Perrenial Glycine

## INTRODUCTION 
Centromeres are essential for chromosome segregation during cell division, yet the underlying proteins and DNA sequences associated with them are remarkably dynamic and rapidly evolving. Centromeric sequences show little conservation even among closely related species as they are prone to rapid mutation and repositioning. This evolutionary plasticity of centromeric sequences contrasts with the consistency of epigenetic regulation and conservation of histone machinery that defines a functional centromere. The dynamic nature of centromere sequences is thought to facilitate cytological processes like centromere repositioning that can drive genomic rearrangements and reproductive isolation during speciation. Increased availability of chromosome-scale and Telomore-2-Telomore genome assemblies for a wide-range of non-model species, fully-haplotype-resolved assemblies, and pangenomes has facilitated deeper study of non-genic regions and genome structral evolution.

In this workflow, I characterized putative centromeric tandem repeats across clades of the Glycine genus. Species with relatively recent common ancestry offer a system for examining the paradoxical combination of rapid sequence evolution and conserved histone function. The Glycine genus is composed of the subgenera Soja and Glycine, which diverged some 10 million years ago (Ma) and differ in life-strategy. Perrenial Glycine represent an extended gene pool for improving annual soybean with traits such as large numbers of seeds per pod, resistance to cyst nematode and fungal pathogens and tolerance to drought and salt stresses (Sherman-Boyles et al., 2014). 

Glycine represents a model system for studying polyploid genome evolution. Diploid Glycine underwent two rounds of whole genome duplication and they share a palaeopolyploidy event with all bean-like (papilionoid) legume species that occurred ~65 Ma with a hare a palaeopolyploidy event with all bean-like (papilionoid) legume species that occurred ~65 Ma. Within the last 350,000 yr there has been a burst of independent allopolyploidy events in the perennial subgenus, with at least eight different allopolyploids (2n = 78, 80) formed from various combinations of eight different diploid genomes.

The genomes for the perrenial Glycine were de novo assembled through a combination of PacBio single molecule real-time sequencing, Illumina sequencing and chromatin conformation capture Hi-C technologies (Methods and Supplementary Fig. 1) and further corrected/improved by integrating previously generated paired bacterial artificial chromosome (BAC) end sequences (BESs) from the same set of accessions.

## RESULTS
To illustrate a proof-of-concept, I took the Glycine max genome assembly, Wm82.gnm6.JFPQ, and ran ULTRA on it to indentify tandem repeats genome-wide.
This produced a JSON file containing data on 1,192,530 tandem DNA repeats! Considering my donwstream workflow wasn't setup to work with JSON files, I opted to convert the ULTRA output to tab-delimited files, using my script [ultra2tsv.v1.sh](https://github.com/CodePipeWrangler/code_jitsu_gh/blob/f6c5753c3ea3fb406f020ca4d49a417e556404b9/genomics/conv_ULTRA-tsv/ultra2tsv.v1.sh). With this script I get TSV files minus the 'sub-repeat' column.

For lightweight data wrangling, I prefer the shell commandline. Shell scripting allows me to rapidly manipulate files, prototype data formats like FASTA or BED files, and validate pipelines without needing to set up full programming environments. To view the distribution of repeat mononer sizes, I can use linux-based shell scripting as follows.

- Histogram of the repeat period distribution (adjust the trailing awk and perl commands respectively for filtering and noise cancellation).
  
    file = ULTRA_FILE_CONVERTED_2_TSV
  
```shell
awk '{print $4}' $file | sort -n | uniq -c | awk '$1>=10 && $2>=60 {print $1 "\t" $2}' | perl -lane 'print $F[1], "\t", "." x int($F[0]/100)'
```

However,for more aesthetic and resolved visualization I can also generate a histogram plot from Python right from commandline as shown below. I simply pass my commands to Python as a string using the '-c' parameter

```shell
    awk '{print $4}' $file | sort -n | uniq -c | sed 's/^ *//' > output.txt
    python -c "import sys; import os; import pandas as pd; import seaborn as sns; import matplotlib.pyplot as plt;
            df=pd.read_csv('output.txt', sep=' ', header=None, names=['Frequency', 'Location'])
            # Create histogram
            plt.figure(figsize=(12, 6))
            plt.bar(df['Location'], df1['Frequency'], color='orange', width=1, edgecolor='black')
            # Labels
            plt.xlabel('ChrORGenXX')
            plt.ylabel('Frequency')
            plt.title('Histogram of Select DNA Repeats Data')
            plt.savefig('plot.png'); 
            print('DONE')"
```
    
On the X and Y axes are repeat or array sizes and their relative frequencies respectively

<figure>
  <img src="https://github.com/user-attachments/assets/adeb46c9-433d-4ce8-9a82-6d57b5553221" alt="Fig.1: glyma.Wm82.gnm6.JFPQ whole genome repeat sizes" width="500" height="300">
  <figcaption>Fig.1 Whole genome distribution of repeat monomer sizes</figcaption>
</figure>

I can also look at distribution of total array sizes to understand more of how the monomers are structured across the genome (adjust the trailing awk and perl commands respectively for filtering and noise cancellation). Here we will create the plot in R.
  
    file = ULTRA_FILE_CONVERTED_2_TSV
  
```shell
awk '{print int($3/1000)}' $file | sort -n | uniq -c | sed 's/^ *//' > glyma.allArrays.txt
```
```R
    library(ggplot2)
    library(dplyr)
    library(tidyr)  # for uncount()

    # Load your aggregated data
    chrom_data <- read.table('glyma.all_Arrays.txt', header = FALSE, sep = "", col.names = c("freq", "array_len"))

    # Expand data by frequency
    expanded_data <- chrom_data %>%
    uncount(freq)

    # Plot histogram + density
    plot <- ggplot(expanded_data, aes(x = array_len)) +
    geom_histogram(aes(y = ..density..), fill = "gold", bins = 50, color = "black") +
    geom_density(color = "firebrick3", size = 1) +
    scale_x_log10() +  # Log transformation
    labs(title = "All Satellite Array Length Distribution",
        x = "Array Length (Kb)",
        y = "Density (log10)")
    print(plot)
    # Save plot
    ggsave("all_arrays_loghistdens_plot.png", plot = plot, width = 10, height = 6, dpi = 300)
```

<figure>
  <img src="https://github.com/user-attachments/assets/8bfd6174-9465-4a5a-a829-ec543997b210" alt="Fig.2: whole genome array sizes" width="500" height="300">
  <figcaption>Fig.2 Whole genome distribution of repeat array sizes</figcaption>
</figure>

On the X and Y axes are array sizes and their relative densities respectively

---

Often a genome contains predomimantly smaller repeats of periods from 1 to ~60 base pairs (bp) that are widely distrubted and show no particular locatization to any chromosomal region (Fig. 1). Repeats larger than 90 bp tend to be distrbuted in hotspots, and thereby potentially representing a centromere. In this case, DNA repeats around 92 bp are highly enriched. The prevailing dogma around centromeric DNA suggest it has several characteristics that can help distinguish it from other types of DNA repeats:

    1) The repeat type is conserved and arranged in large arrays on the chromosomes
    2) Localization of a repeat type to discrete genetic regions on many or all chromosomes
    3) High abundance of the repeat in discrete regions.
    4) Low diversity of repeat type in discrete regions
    5) The repeat type is larger than 60 bp 

---

#### I created some ridgeline plots in R to visualize the whole genome, tandem repeat distribution. 

I began by plotting arrays of all repeat periods (i.e. monomer sizes) by chromosome. To do this I used the follwoing shell commands to create the input file for the histrogram.

```shell
for i in {01..20}; do echo Gm$i ; awk ''/'Gm'$i'/ {print int($2/1000000)}' $file | sort -n | uniq -c | awk '{print $0" "'$i'}' | sed 's/^ *//' >> test.txt; done
```

Then for selected repeat sizes of 91, 92, and 91-92 bp, since I know these represent the Glycine max and Glycine soja centromeric repeat sizes

```shell
for i in {01..20}; do echo Gm$i ; awk ''/'Gm'$i'/ && $4==91 {print int($2/1000000)}' $file | sort -n | uniq -c | awk '{print $0" "'$i'}' | sed 's/^ *//' >> test2.txt; done
```
```shell
for i in {01..20}; do echo Gm$i ; awk ''/'Gm'$i'/ && $4==92 {print int($2/1000000)}' $file | sort -n | uniq -c | awk '{print $0" "'$i'}' | sed 's/^ *//' >> test3.txt; done
```

```shell
for i in {01..20}; do echo Gm$i ; awk ''/'Gm'$i'/ && $4>=91 && $4<=92 {print int($2/1000000)}' $file | sort -n | uniq -c | awk '{print $0" "'$i'}' | sed 's/^ *//' >> test4.txt; done
```

*Next I'll slide over into my R console and plot...*

```R
#Import libraries
library(ggplot2)
library(ggridges)
library(viridis)

# Move to working directory
setwd("/PATH/TO/YOUR/FILES")

# Load the data
chrom_data <- read.table('test.txt', header = FALSE, sep = " ", col.names = c("frequency", "position", "chromosome")) # the 3rd column will change according to the data. It is used for distinguishing between plots.
chrom_data$chromosome <- factor(chrom_data$chromosome)

# Plotting
# Create a ridgeline plot using ggridges
plot <- ggplot(chrom_data, aes(x = position, y = chromosome, height = frequency, fill = after_stat(x), group = chromosome)) +
  geom_density_ridges_gradient(stat = "identity", scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Frequency", option = "C") +
  theme_minimal() +
  labs(title = "Glycine max (Wm82.gnm6) All DNA Repeats",
       x = "Genomic Position",
       y = "Chromosome")
print(plot)
ggsave("Glycine_max_chromosome_ridgeline_plot.png", plot = plot, width = 10, height = 6, dpi = 300)

```
    
![Fig.3: Glycine_max_chr_ridgeline_plot](https://github.com/user-attachments/assets/71750bab-614c-4821-b520-6dcc19e787ff)

Below are ridgeline plots from selected chromosomes to illustrate that 91 and 92 bp satellites are concentrated in single hotspots across chromosomes.

<figure>
  <img src="https://github.com/user-attachments/assets/6bef897d-bfe3-471e-99ea-8efc24258dbd" alt="Fig.4: Gm2 centromeric repeat dist." width="300" height="300">
  <img src="https://github.com/user-attachments/assets/6962430f-b41b-465a-878f-a14fe12c0b9b" alt="Fig.5: Gm3 centromeric repeat dist." width="300" height="300">
  <img src="https://github.com/user-attachments/assets/1876e37c-7063-47c1-bbe1-08dbae4b54df" alt="Fig.6: Gm5 centromeric repeat dist." width="300" height="300">
  <img src="https://github.com/user-attachments/assets/a6052064-92cd-42d5-a8eb-60ff9791d490" alt="Fig.7: Gm9 centromeric repeat dist." width="300" height="300">
  <img src="https://github.com/user-attachments/assets/2a8a60cf-9f0e-409f-bd28-58964edb31f7" alt="Fig.8: Gm13 centromeric repeat dist." width="300" height="300">
  <img src="https://github.com/user-attachments/assets/65b48cec-5dd3-400c-a875-ddb9ab32efc1" alt="Fig.9: Gm20 centromeric repeat dist." width="300" height="300">
</figure>

We can see that the 91-92 bp repeats satisfy several of the criteria for a centromeric repeat. Many of the largest arrays in Wm82.a6 correspond to these repeats, and in this case there were no consistent hotspots of other large arrays to compete for the role of centromeric repeat. 

<figure>
  <img src="https://github.com/user-attachments/assets/615e135c-dec2-4627-8245-e1e75eda486b" alt="Fig.10: 91-bp whole genome array sizes" width="300" height="300">
  <img src="https://github.com/user-attachments/assets/02cf7bea-6b36-4578-8018-3f09bffd5c02" alt="Fig.11: 92-bp whole genome array sizes" width="300" height="300">
</figure>

I can also take data from ULTRA and run the code below to look at the distribution of 91-92 bp repeats along chromosomal tracks
Using the plotting file (i.e. test4.txt) and methodolgy I showed above for ridgeling plots ...

![Fig 12: Glycine_max_91-92bp_chromosome_ridgeline_plot](https://github.com/user-attachments/assets/259253fa-4266-4a90-bb0e-dd571f48825e)

Furthermore, previous research on this subject gives support that these 91 and 92 bp represent CentGm-1 and 2, the two types of *Glycine max* centromeric repeats.

---
There is another program, RepeatObserver, to assist with characterization of tandem DNA repeats. The creators of *RepeatObserver* found that centromeres consistently occur in regions that have the least diversity in repeat types (i.e. one or a few repeated sequences are present in very high numbers). Taking advantage of this information, a genomic Shannon diversity index is used to predict centromere locations in several other chromosome-scale genome assemblies. The Fourier spectra produced by RepeatOBserver can help visualize historic centromere positions, potential neocentromeres, retrotransposon clusters and gene copy variation. Identification of patterns of split and inverted tandem repeats at inversion boundaries suggests that at least some chromosomal inversions or misassemblies can be predicted with RepeatOBserver.

RepeatObserver outputs many files that can be helpful for characterizing centromeres. I'll show a few here that were particularly helpful for me*. RepeatObserver uses several metrics to predict the centromeric content and boundaries for you. In a model system like soybean, these plots often show high-resolution.

RepeatObserver: Centromere Histograms, Shannon diversity plots, and centromere prediction plots

<img src="https://github.com/user-attachments/assets/98d50066-3643-4d14-ad6b-b43b37fa484f" alt="Fig.13: RepeatObserver: Centromere Histograms" width="300" height="300">
<img src="https://github.com/user-attachments/assets/7d0afc41-d036-4618-8cfe-5d26d05ae0d2" alt="Fig.14: RepeatObserver: Shannon diversity plots" width="300" height="300">
<img src="https://github.com/user-attachments/assets/b8b5e242-e5b9-4c6f-8b2d-7f07ae89646f" alt="Fig.15: RepeatObserver: Centromeric predictor comparisons" width="300" height="300">

Let's take a look at a few chromosomes in particular. Here is a Fourier transform plot from RepeatObserver for Chr.2. Such plots in RepeatOBserver provide a visualization of a chromosome's DNA sequence in the frequency domain, generated by applying a Fourier transform to a "DNA walk", whis is a simplified representation of the AT/CG patterns. This visualization facilitates identification of regions with repetitive patterns, particularly centromeres, by looking for peaks in the plot that represent the dominant repeat lengths present at that location on the chromosome. 

<img src="https://github.com/user-attachments/assets/9ef00dbb-2f60-45be-b9ff-355571135417" alt="Fig.16: Gmax_H0-AT_Chr2_All_spec_bp35_2000seq1_10959TRUE" width="300" height="300">
<img src="https://github.com/user-attachments/assets/e1b78e0f-5145-4b84-9b51-179d96c668df" alt="Fig.17: Gmax Chr2 91-92 bp repeats" width="300" height="300">

We can see that only one type of repeat is dominant at the peak on Chr.2. From ULTRA, it is evident 91-92 bp repeats are enriched in the genome and localized to specific regions on each chromosome. This is characteristic of a centromeric repeat. Looking at the distribution of 91 and 92 bp repeats on this chromosome, it lines up directly with the peak in the fourier transform plot. 

---

Then I wanted to extract the putative centromeric repeats and analyze the sequence further. Are the repeats alike? transposons? How mych do they vary? from location to location in the genome? on the same chromosome? Many questions can come to mind, so I grabbed all of the results for putative centromeric repeats and did phylogenetic analyses to get some answers. For *Glycine max* these repeats are similiar enough from place to place that we can pool them all together and still proceed with the analysis without issue. However, there are many cases where it benefits clarity in the results to take repeats from bins of interest and then analyze those. For instance, I may pick 1 Mbp around the peak of 91 bp repeats and pool sequences from every chromosome to determine if they are conserved in the genome. Although centromeric sequences can be highly variable between species, it is expected that variation is less than 5-10% within a species. This is evident in *Glycine max*, where the 91 and 92 bp repeats from Wm82, Lee, and FiskebyIII are nearly identical.
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

Sometimes you may see a repeat type on a chromosome or at a neo-centromere that looks to be significantly different in composition from others in the set. It could be there are more than one repeat type in the data, even possibly of the same size. However, chromosomal restructuring throughout evolution can invert DNA sequences or you could be looking at the other strand. To investigate this we can take the reverse complements of the consensus sequences and check their similarity. In this case, the reverse complements are 91% and 92% similar to the Wm82.gnm6.cons-1 and 2 sequences respectively.

---

* INSERT blast results for both repeat types: where do the two types of repeats localize to? same chromosomes? same copy number?*


*two types of repeats in Glycine max.*

CentGm-1 and 2 show large variation in copy number. On most chromosomes CentGm-1 is dominant, while others show CentGm-2 or codominance.



