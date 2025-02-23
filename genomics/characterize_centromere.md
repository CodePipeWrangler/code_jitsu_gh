# Objective: characterize centromeric repeats in Glycine


## INTRODUCTION 
While centromeres are essential for chromosome segregation during cell division, the underlying DNA sequences associated with centromeres are remarkably dynamic and rapidly evolving. Centromeric sequences show little conservation even among closely related species as they are prone to rapid mutation and repositioning. This evolutionary plasticity of centromeric sequences contrasts with the consistency of epigenetic regulation and conservation of histone machinery that defines a functional centromere. The dynamic nature of centromere sequences is thought to facilitate cytological processes like centromere repositioning that can drive genomic rearrangements and reproductive isolation during speciation.

Increased availability of chromosome-scale and Telomore-2-Telomore genome assemblies for a wide-range of non-model species, fully-haplotype-resolved assemblies, and pangenomes has facilitated deeper study of non-genic regions 
and genome structral evolution.

In this workflow, I characterized putative centromeric tandem repeats across clades of the Glycine genus. Species, with relatively recent common ancestry, offer a system for examining the paradoxical combination of rapid sequence evolution and conserved histone function. The Glycine genus is composed of the subgenera Soja and Glycine, which diverged some 10 million years ago (Ma) and differ in life-strategy. Perrenial Glycine represent an extended gene pool for imporving annual soybean with traits such as large numbers of seeds per pod, resistance to cyst nematode and fungal pathogens and tolerance to drought and salt stresses (Sherman-Boyles et al., 2014). Glycine represents a model system for studying polyploid genome evoltuion. dip;oid Glycine undrwent two rounds of whole genome duplication and they share a palaeopolyploidy event with all bean-like (papilionoid) legume species that occurred ~65 Ma with a hare a palaeopolyploidy event with all bean-like (papilionoid) legume species that occurred ~65 Ma. Within the last 350,000 yr there has been a burst of independent allopolyploidy events in the perennial subgenus, with at least eight different allopolyploids (2n = 78, 80) formed from various combinations of eight different diploid genomes.

The Glycine genomes for the perrenial Glycine were de novo assembled through a combination of PacBio single molecule real-time sequencing, Illumina sequencing and chromatin conformation capture Hi-C technologies (Methods and Supplementary Fig. 1) and further corrected/improved by integrating our previously generated paired bacterial artificial chromosome (BAC) end sequences (BESs) from the same set of accessions.

## RESULT
To illustrate a proof-of-concept, I took the Glycine max genome assembly, Wm82.gnm6.JFPQ, and ran ULTRA to indentify tandem repeats genome-wide.
This produced a JSON filke containing data on 1,192,530 tandem DNA repeats! Considerign my donwstream workflow wasn't comfortable with JSON files, I opted to convert the ULTRA output to tab-delimited files, using my script *ultra2tsv.sh*. With this script I get TSV files minus the 'sub-repeat' data, which I do not find of benefit to this research in partciluar.

To view the distribution of repeat mononer sizes, I can use linux-based commandline tools as follows:

### Lets see a histogram of the repeat period distribution (adjust the trailing awk and perl commands respectively for filtering and noise cancellation). On the left column are repeat sizes followed by their relative frequency to toehr repeats. Often your genome contains mostly
smaller repeatrs of periods from 1 to ~60 that are widely distrubted and show no partuclar locatization to any chromosomal region. Repeats lartger than 90 bp tend to be distrbuted in hotspots, potentially representing a centromere.
    awk '{print $4}' $file | sort -n | uniq -c | awk '$1>=10 && $2>=60 {print $1 "\t" $2}' | perl -lane 'print $F[1], "\t", "." x int($F[0]/100)'
----
![Freq. of Wm82.v6 Tandem DNA Repeats by monomer size](/Users/bjordan/Desktop/Screenshot 2025-02-23 at 2.35.40 PM.png)


The largest tandem repeat arrays were enriched for satellite sequences of particular periods that were conserved across chromosomes within each species


