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

Taking a quick look using the Linux commands below, a histogram of the repeat period distribution is generated (adjust the trailing awk and perl commands respectively for filtering and noise cancellation). 
    
    awk '{print $4}' $file | sort -n | uniq -c | awk '$1>=10 && $2>=60 {print $1 "\t" $2}' | perl -lane 'print $F[1], "\t", "." x int($F[0]/100)'
    
On the left column are repeat sizes followed by their relative frequency to other repeats:

![glyma Wm82 gnm6 JFPQ large](https://github.com/user-attachments/assets/adeb46c9-433d-4ce8-9a82-6d57b5553221)


Often your genome contains mostly smaller repeats of periods from 1 to ~60 that are widely distrubted and show no particular locatization to any chromosomal region. Repeats larger than 90 bp tend to be distrbuted in hotspots, potentially representing a centromere. I can see that DNA repeats around 92 bp are highly enriched in this dataset. The prevailing ideas in plant genetics literature are that centromeric DNA has several characteristics that can help distinguish it from other types of DNA repeats:
    1) The repeat type is conserved and arranged in large arrays on the chromosomes
    2) Localization of a repeat type to discrete genetic regions on many or all chromosomes
    3) High abundance of the repeat in discrete regions.
    4) Low diversity of repeat type in discrete regions
    5) The repeat type is larger than 60 bp 

![glyma Wm82 gnm6 JFPQ Allarrays](https://github.com/user-attachments/assets/85773ba2-39cd-4798-86d6-178b3518367f)
We can see that the 91-92 bp repeats satisfy several of these criteria. 1) Many of the large arrays in Wm82.a6 correspond to these repeats. Although, not all large arrays are made up of centromeric repeats. Previous work on this subject gives confidence that these repeats of periods 91-92 represent CentGm-1 and 2, being the two types of *Glycine max* centromeric repeats.
![glyma Wm82 gnm6 JFPQ 91-92arrays](https://github.com/user-attachments/assets/dc6062a2-4801-4dab-a190-820e973b2190)
![glyma Wm82 gnm6 JFPQ Altarrays](https://github.com/user-attachments/assets/d943cd18-fd32-4b87-ba30-d8ace1ac7e74)


I can also take data from ULTRA and run the loop below to look at the distribution of 91-92 bp repeats along a gene track:
    file=ultra.glyma.Wm82.gnm6.JFPQ.p1000.tsv

    for i in {01..20}; do echo Gm$i ; 
        ref=`echo $file | perl -pe 's/ultra\.(.+)\.p1000\.tsv/$1/'`                                    
        awk ''/'.Gm'$i'/ && $4==91 {print int($2/1000000)}' $file |  #note 'Gm' is a specific chromosomal identifier for this species. 
        sort -n | uniq -c | sed 's/^ *//' > output$i.txt          
        locs=`cat output$i.txt | wc -l`                                                              
        export var11=$ref         
        export VAR31=$i 
        python -c "import sys; import os; import pandas as pd; import seaborn as sns; import matplotlib.pyplot as plt;
            var11 = os.getenv('VAR11');var21=os.getenv('VAR21');var31=os.getenv('VAR31')
            df1=pd.read_csv(f'output{var31}.txt', sep=' ', header=None, names=['Frequency', 'Location'])
            # Create histogram
            plt.figure(figsize=(12, 6))
            plt.bar(df1['Location'], df1['Frequency'], color='orange', width=1, edgecolor='black')
            # Labels
            plt.xlabel(f'Gm{var31}')
            plt.ylabel('Frequency')
            plt.title('Histogram of 91-bp Repeats')
            plt.savefig(f'{var11}.91.{var31}.png'); print('DONE')"
    done

Before moving on, let's revisit the grpah of epriod sizes above. More support is needed to confidently say 91-92 bp repeats represent the centromeric repeat. What if there are other, less obvious repeats that underly the centromeric machinery. Hence, without prior knowledge of the centromeric repeat type (e.g. its length), how else would we decipher which is a centromeric repeat?

To investigate this I ran another program, RepeatObserver, to assist with further characterization of tandem DNA repeats. The authors of *RepeatObserver* found that centromeres consistently occur in regions that have the least diversity in repeat types (i.e. one or a few repeated sequences are present in very high numbers). Taking advantage of this information, a genomic Shannon diversity index is used to predict centromere locations in several other chromosome-scale genome assemblies. The Fourier spectra produced by RepeatOBserver can help visualize historic centromere positions, potential neocentromeres, retrotransposon clusters and gene copy variation. Identification of patterns of split and inverted tandem repeats at inversion boundaries suggests that at least some chromosomal inversions or misassemblies can be predicted with RepeatOBserver.

RepeatObserver outputs many files that can be helpful for characterizing centromeres. I'll show a few here that really helped distinguish where the centromere is located.

RepeatObserver uses several metrics to predict the centromeric content and boundaries for you. In a model system like soybean, these plots often show high-resolution.

![Gmax_H0-AT_Histograms](https://github.com/user-attachments/assets/98d50066-3643-4d14-ad6b-b43b37fa484f)
![Gmax_H0-AT_Shannon_div](https://github.com/user-attachments/assets/7d0afc41-d036-4618-8cfe-5d26d05ae0d2)
![Gmax_H0-AT_centromere_prediction_comparisons](https://github.com/user-attachments/assets/b8b5e242-e5b9-4c6f-8b2d-7f07ae89646f)

*cluster reaults and alignment of many 91-92 bp sequences.


### Extract data for specific localized repeat sizes and Create fasta file
    x='156'; pos1=48000000 ; pos2=52000000; i=1
    for file in *tsv; do
    for sub in `echo $x`; do
    awk -v CUT="$sub" -v START="$pos1" -v STOP="$pos2" '$4==CUT && $2>=START && $2<=STOP && '/'chr'$i'/ {print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $file >> B97.155-156.chr1.centCANDI-ONLY.fn ;
    done; 
    done
    
*alignments of putative and neo-centromeric 91-92 Bp repeat monomers.
    - 91 alignments 
![glyma gnm6 91 ALL aln](https://github.com/user-attachments/assets/638439fd-60dd-4550-863a-a4346b6ea0b8)
Short alignments of tandem repeat monomers can be misleading, given ULTRA defines what is a repeat monomer somewhat arbitrarily. To assess whether real differencs are beign seen in the alignment, I increased the size of every repeat by four times (see code below) and realligned them.

     awk '/^>/ {h=$0; print h} !/^>/ {for(i=0; i<4; i++) print $0}' "$file" > "${file%.fn}.4X.fn"

<img width="804" alt="glyma gnm6 91 ALL 4X aln" src="https://github.com/user-attachments/assets/d14fbb7d-1f2a-499b-b297-587a3e3a5428" />

The same phenomenom was evident from aligning 92 bp repeats from ULTRA. Trimming the 4X alignment visually to one unit of the repeat, I then made consensus sequen ces for each using the tool [EMBOSS Cons](https://www.ebi.ac.uk/jdispatcher/msa/emboss_cons?stype=protein) .

Consensus sequences:
>glyma.gnm6.JFPQ.cent_91.cons
AAAAAGTTATTGTCGTTTAAATTTGCTCAGAGCTTCATTTTTCAATTTCGAGCGTCTCGATATATTACGGGACTCAATCAGACATCCAATT
>glyma.gnm6.JFPQ.cent_92.cons
AAAAGTTATGACCATTTGAATTTCTCGAGAGCTTCCGTTGTTCAATTTCGAGCGTCTCGATATATTATGCGCCTGAATCGGACATCCGAGTG

Sometimes you may see a repeat type on a chromosome or at a neo-centromere that looks to be signiifcantly different in compostion from others in the pool. It could be there are more than one repeat type in the data, even possibly of the same size. However, chromosomal restructruing throughout evolution can invert DNA sequences or you could be looking at the other strand. To investigate this we can take the reverse complements of the consensus sequences and check their similarity.
   
>glyma.gnm6.WGD.cent_91.30percent.cons
AATTGGATGTCTGATTGAGTCCCGTAATATATCGAGACGCTCGAAATTGAAAAATGAAGCTCTGAGCAAATTTAAACGACAATAACTTTTT
>glyma.gnm6.WGD.cent_92.30percent.cons
CACTCGGATGTCCGATTCAGGCGCATAATATATCGAGACGCTCGAAATTGAACAACGGAAGCTCTCGAGAAATTCAAATGGTCATAACTTTT


Alignment statistics for match glyma.gnm6.JFPQ.cent_91.cons vs. reverse complement
Score	Expect	Identities	Gaps	Strand
165 bits(182)	2e-46	91/91(100%)	0/91(0%)	Plus/Minus
Query  1   AAAAAGTTATTGTCGTTTAAATTTGCTCAGAGCTTCATTTTTCAATTTCGAGCGTCTCGA  60
           ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  91  AAAAAGTTATTGTCGTTTAAATTTGCTCAGAGCTTCATTTTTCAATTTCGAGCGTCTCGA  32

Query  61  TATATTACGGGACTCAATCAGACATCCAATT  91
           |||||||||||||||||||||||||||||||
Sbjct  31  TATATTACGGGACTCAATCAGACATCCAATT  1

Alignment statistics for match glyma.gnm6.JFPQ.cent_92.cons vs. reverse complement
Score	Expect	Identities	Gaps	Strand
167 bits(184)	6e-47	92/92(100%)	0/92(0%)	Plus/Minus
Query  1   AAAAGTTATGACCATTTGAATTTCTCGAGAGCTTCCGTTGTTCAATTTCGAGCGTCTCGA  60
           ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  92  AAAAGTTATGACCATTTGAATTTCTCGAGAGCTTCCGTTGTTCAATTTCGAGCGTCTCGA  33

Query  61  TATATTATGCGCCTGAATCGGACATCCGAGTG  92
           ||||||||||||||||||||||||||||||||
Sbjct  32  TATATTATGCGCCTGAATCGGACATCCGAGTG  1

In this case, the reverse complement is highly similar to the CentGm-1 and 2 sequences.

*blast results for both: where do the two types of repeats localize to? same chromosomes? same numbers?*

*two types of repeats in Glycine max.*




