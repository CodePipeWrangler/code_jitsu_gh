# Shell Jutsu: The Way of the Command Line

Just some of my command-line jutsu for data wrangling and analysis. Regex throws, awk strikes, and tarball flips—perfect for wrangling massive datasets with minimal motion.

*"Absorb what is useful, discard what is not, add what is uniquely your own." - Bruce Lee*

### File and directory management
ls PATH/TO/DIR – List files

cd PATH/TO/DIR – Change directory

basename "$PWD" – Get current directory name only

mkdir foo && cd "$_" – Make directory and move into it

mkdir -p /foo/bar && cp myfile.txt $_ – Make directory and copy files into it

rsync -avu --delete "A/" "B/" – Sync folders and delete or update files

### Regex utilities
These are some Regex tools that really help me with my work. There are many times when I need to match a specific pattern to manipulate a file or extraxct information from it.
#### Match any SET of alphabetic or numeric characters and also the '_' (underscore). Particular special characters can be added to this schema.
    [A-Za-z0-9_]+
Match either case where regex ends in 'fn' or 'fan'
- fna?
    - Matches fn followed by a (for fna) OR fn without a (for fa)

List files in a diretory
- $ ls PATH/TO/DIR	

Change the working directory
- $ cd PATH/TO/DIR	

Get the current directory without the full PATH
- $ basename "$PWD"

See statistics of currently running processes
- $ top

Transpose a column to a row
[FEED] |  tr '\n' ' '

Tar zip a file(s)
- $ tar -czf  name-of-archive.tar.gz /path/to/directory-or-file /path/to/another-directory-or-file …

See contents of zipped file without extraction
- $ tar -tf filename.tar.gz
- $ vim filename.tar.gz
    - This method allows the user to also see the file contents.
- $ less filename.tar.gz
- $ zcat filename.tar.gz
    - To grep something - $ zgrep test  filename.tar.gz
    - To check difference between files - $  diff filename1.tar.gz  filename2.tar.gz


Basic arithmetic including variables
- $ echo $((1+1))
    - 2
- some_var=2
- $ my_var=$((1+1+some_var)) ; echo my_var
    - 4
- 
Access the manual page for a command
- $ man cd
 
word count: print the number of lines
- $ wc -l

Execute command from history
- see "history expansion" LESS='+/^HISTORY EXPANSION' man bash
- In bash, just 
	$ !636 
- run last executed command: 
	$ !!
- run two commands back: 
	$ !-2
- repeat the last argument of the last command (this opens the saved file)
	$ echo this text goes into a file > /tmp/afile.txt
	$ cat !$

Mkdir and cd to it
mkdir foo && cd "$_"
- $_is a special parameter that holds the last argument of the previous command. The quote around $_ make sure it works even if the folder name contains spaces.
- Why use double quotes?
	In some shells, such as zsh, the 		double quotes surrounding the $_ 	are not necessary even when the directory name contains spaces. They are required for this command to work in bash, however.

Mkdir and cp files to it (one liner):
- $ mkdir -p /foo/bar && cp myfile.txt $_

Print a new line:
- $ echo ""
- $ echo $'hello\nworld'
	prints
		hello
		world
- $'' strings use ANSI C Quoting:
	Words of the form $'string' are treated specially. The word expands to string, with backslash-escaped characters replaced as specified by the ANSI C standard.

Print two (or more?) blank lines
- $ echo -e '\n'
- Add more \n to increase the number of blank lines OR…
Print a FEW blank lines (using bash or ksh93 or zsh)
- $ printf '%.0s\n' {1..3}

Show only current directory name not full path
- $ ${PWD/*\//}
- $ ${VAR/pattern_to_find/pattern_to_replace}`

Rename files by subbing string
- $ rename  [pattern] [replacement] file

Chain linux commands executing subsequent commands only if previous return with zero exit status:
- $ Command1 && command2

Run loop on specified files:
- $ for file in *; do; done

Get lines matching a string
- grep 'STRING' filename 

Get line matching a string an X lines around it using grep
- Grep has the following options. See the man page for more information:
- -A [num] Print num lines of trailing context after each match. See also the -B and -C options.
- -B [num] Print num lines of leading context before each match. See also the -A and -C options.
- -C [num] Print num lines of leading and trailing context surrounding each match. The default is 2 and is equivalent to -A 2 -B 2. Note: no whitespace may be given between the option and its argument.
- 
Get lines matching one or more patterns and send to a file:
- $ grep ' P \| CA ' file > new_file
- $ grep -E ' (P|CA) ' file > new_file
- $ awk '/ P / || / CA /' file
- On Mac OS Ventura:
    - $ grep -e ' CA ' -e ' P ' all.pdb > CA.pdb 
    - From the man page of grep: -e pattern, --regexp=pattern ; Specify a pattern used during the search of the input: an input line is selected if it matches any of the specified patterns. This option is most useful when multiple -e options are used to specify multiple patterns, or when a pattern begins with a dash (‘-’).

Search files in dir for pattern and save to file:
- $ grep 'P|CA' -Er /path_to_your_dir/ > /tmp/grep.log
- If you need case insensitive, replace -Er to -Eri. In file /tmp/grep.log you will see path to file and matched string. 
- if you need search in files with specific extension then write something like:
    - $ grep 'P|CA' -Er --include=*.php /path_to_your_dir/ > /tmp/grep.log

Make file executable 
$ chmod +x file.sh

Get number of commandline arguments
$ $#

Redirect std-error and std-out to a file
- [command] 2>&1

- > redirects output to a file, overwriting the file.
- >> redirects output to a file appending the redirected output at the end.
- Standard output is represented in bash with number 1 and standard error is represented with number 2. They are separate, so the user can redirect them to different files.
- 2>&1 redirects the standard error to the standard output so they appear together and can be jointly redirected to a file. (Writing just 2>1 would redirect the standard error to a file called "1", not to standard output.)
More, less, and most
- See https://askubuntu.com/questions/1191862/what-is-the-difference-between-more-and-less-commands
Count how many columns in file:
- $ cat file | awk '{print NF}'

Remove a column using awk
- $ awk '!($10="") {print}' FILE

Search file for pattern in particular col and print
- $ awk '$1~/Wm82.gnm2.ann1/ {print}'

Transpose tab-sep line to column
- $ tr '\t' '\n' file

Get max number in column:
- $ cat file.tsv | awk -v max=0 '{if($1>max){want=$2; max=$1}}END{print want}';
- This works on two column data where you 'want' the maximum value over all other values. 'Max' is a placeholder for the current max value before analyzing the entire dataset.

Get the nearest rounded number
- $ var=2.5
- $ echo $var | awk '{print int($1+0.5)}'

*************************************************************************************
******************************STOPPED_HERE_EDITING*****************************
*************************************************************************************

Use shell variable in awk regex pattern
This must be accomplished in a different manner than passing an awk variable to print. I presented this concept below encased in a for loop that runs an awk process for 20 chromosomes denoted as Gm01, etc.

for x in {01..20}; do echo Gm$x ; awk 'FNR>5 && '/'.Gm'$x'/ && $3>=90 {print $9}' out.blastn.Gmr12_rep.glyma.*ISU* | histogram -n | grep -v '0' ; echo; echo; done

Case-insensitive regex search for negative of pattern at start of pattern
/^[^PATTERN]/I

Get average of a column with Awk:
$ head -15 *Wm82_I*gnm4*.gff | awk '{print $9}'| grep -o -P '(?<=value=).*(?=;ID)'
- So from the column below: 
* value=191;ID=191
* value=146;ID=337
* value=70;ID=407
* value=21;ID=428
* value=215;ID=643
* value=65;ID=708
* value=123;ID=831
* value=456;ID=1287
* value=23;ID=1310
* value=51;ID=1361
* value=174;ID=1535
* value=110;ID=1645
* value=173;ID=1818
* value=177;ID=1995
* value=199;ID=2194
	I get only the actual values designated by value.
- Now I pipe that to… | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }'
- Calculate the median with… | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'
- OR use median.awk program

Get average and standard deviation (population) of a column with Awk
$ awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} 
          END {for (i=1;i<=NF;i++) {
          printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
         }' file.dat >> aver-std.dat

Get standard deviation (sample) of a column with Awk
$ awk '{sum+=$0;a[NR]=$0}END{for(i in a)y+=(a[i]-(sum/NR))^2;print sqrt(y/(NR-1))}' $file

	- "The two forms of standard deviation are relevant to two different types of variability. One is the variability of values within a set of numbers and one is an estimate of the variability of a population from which a sample of numbers has been drawn.
The population standard deviation is relevant where the numbers that you have in hand are the entire population, and the sample standard deviation is relevant where the numbers are a sample of a much larger population.
For any given set of numbers the sample standard deviation is larger than the population standard deviation because there is extra uncertainty involved: the uncertainty that results from sampling. See this for a bit more information: Intuitive explanation for dividing by
𝑛−1
when calculating standard deviation?
For an example, the population standard deviation of 1,2,3,4,5 is about 1.41 and the sample standard deviation is about 1.58."
                            - from https://stats.stackexchange.com/questions/485326/confused-when-to-use-population-vs-sample-standard-deviation-in-engineering-test

Cut columns except:
$ cut -d$'\t' -f 1-10

Extract parts of a string and create a table of the data
This takes the stated string and parses it to look like the output.

% echo 'High homology BLAST matches to the glyd3.G1403.gnm1.Chr11_27836031_1661_454 nucletide periodic sequence' && echo '-----------------------------------------------------------------------------------------------------------' && 
echo 'spec.\tchr.\tstart.pos.\tlen.\tperiod.\tperc.iden\tmatch.len' && awk -v OFS='\t' '$1~/glyd3.G1403.gnm1.Chr11_27836031_1661_454/ && $3>=80 {print $2,$3,$4}' *blast*txt | perl -pe 's/(\w+)\.\w+\.\w+\.(Chr\d\d)\_(\d+)\_(\d+)\_(\d+)/$1\t$2\t$3\t$4\t$5\t/' | sort -k1 -k2 

Convert json file to tabular format with preset columns
for file in gly*json; do echo $file; cat $file | perl -pe 's/"//g; s/,//; s/\{//; s/\}//; s/://' | awk -v ORS="" '$1~/SequenceName/ {print "\n"} $1~/Start|Length|Period|Score|Substitutions|Insertions|Deletions|Consensus|Sequence/ {print $2 "\t"}' > $file.tsv; echo; done 

Update a directory from another
rsync -avu --delete "/home/user/A/" "/home/user/B"
* -a Do the sync preserving all filesystem attributes
* -v run verbosely
* -u only copy files with a newer modification time (or size difference if the times are equal)
* --delete delete the files in target folder that do not exist in the source

Copy multiple files
cp /home/usr/dir/{file1,file2,file3,file4} /home/usr/destination/
- The syntax uses the cp command followed by the path to the directory the desired files are located in with all the files you wish to copy wrapped in brackets and separated by commas.
- Make sure to note that there are no spaces between the files. The last part of the command, /home/usr/destination/, is the directory you wish to copy the files into.
- or if the all the files have the same prefix but different endings you could do something like this:
cp /home/usr/dir/file{1..4} ./
- Where file1,file2,file3 and file4 would be copied.
- 


Install program on Ceres via miniconda while creating env
conda create -n mummer4 -c conda-forge -c bioconda perl-bioperl-core mummer4 gnuplot

Add prefix to filenames
for f in * ; do mv -- "$f" "PRE_$f" ; done
	ex.
		for f in * ; do mv -- "$f" "phavu.$f" ; done

Remove suffix from file names
for f in *; do mv "${f}" "${f/_SUF/}"; done
	ex.
		for file in *.tmp; do mv "${file}" "${file/.tmp/}"; done

Vim indent/unindent a selection
way is to select a block and insert an indent at the beginning of the line using this sequence:
1. ctrl+V + arrow keys to select the block.
2. I to switch to insert mode such that the inserted text is inserted at the beginning of the selection in each line in the selected block.
3. ctrl+T to increase the indent or ctrl+D to decrease the indent. You can add any number of indents like this. Note: The indentation will be seen only the first line of the block, but when insert mode is exited the indentation will be replicated on all the lines in the block.

Vim side by side / diff views
- See https://unix.stackexchange.com/questions/1386/comparing-two-files-in-vim
- Vimdiff file.txt file2.txt
- Split screen mode: vim -0 file1 [file2 …]
    - Then rune : diffthis OR FOR OFF diff off
- OR : vs otherfile
- ctrl+w h or l to swap screens
- : Diffthis ti turn on diff mode in either screen
- Vim -d file1 [file2…] is equal to calling vimdiff directly
- :vs otherfile literally compares two files in vim, not vimdiff
Vim remove last search pattern
:nohlsearch OR :not

Vim search and replace within selection
Highlight text then enter…
:%s/\%VSEARCH/REPLACE/g
%s/\%Vmod0/mod2/g

***Read list items into  loop from text file***:
for file in `cat test.txt`; do echo $file; done

Linux: check existence of dir and variable
[-d "path/to/dir']
[-z ${var+x} ] # if variable does not exist, then null. Otherwise the variable is subbed by x.

Notes on scripting
a shebang MUST be present (e.g., #!/bin/sh, if POSIX compliant)

Get GitHub CLI version installed
- $ gh version

Sort by one column and then the next
sort -k 1,1 -k2,2n file

Remove duplicates from a fasta file
Ml seqkit
seqkit rmdup -s < in.fa > out.fa

Remove line breaks from sequences in a fasta file
$ awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' input.fasta

Find files created after date and with a suffix (e.g. '.json')
$ find . -type f -name "*.json" -newermt "2024-08-01" -exec ls -l {} +

Find filed created within a certain date range and matching a pattern
$ find . -type f -name "*.json" -newermt "2024-09-17" ! -newermt "2024-09-19" -exec ls {} +  

### Find files except that matching a PATTERN
find . -maxdepth 1 ! -name "*fna" -exec rm -r {} + # The option maxdepth restricts the command to only run in the current directory

Find multiple file patterns at once (could this be simpler?)
find . \( -name "G*glys*gff" -o -name "G*glyd*gff" \)
Grep 'PATTERN|PATTERN'
- Note grep regex is different,like perl, than bash regex
- 
Add a folder or file to your PATH variable
- vim ~/.bashrc
Add this line to the file: 
- export PATH="/path/to/your/folder:PATH"
Save the file. Then call the following line or restart the session to apply the changes to your current session.
- Source ~/.bashrc

Split a string into chunks of 6 and append a prefix to it
$ echo Hello World | fold -w6 | sed -e 's/^/chunk_/'

Use desktop calculator command
$ dc
https://www.computerhope.com/unix/udc.htm#:~:text=If%20called%20from%20the%20top,command%20causes%20dc%20to%20exit.

Show how significant the average lengths of sequences in a Fasta file (e.g. set of chromosomes) are using seqlen.awk from my bin and awk programming
$ for file in gly*/*main.fna; do echo $file; seqlen.awk $file | cut -f2 | awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}}                     
          END {for (i=1;i<=NF;i++) {
          printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
         }' > data.txt ; cat data.txt; 
        while IFS=" " read -r ref1 ref2 remainder; do
                high=`echo $ref1+$ref2 | bc` ; low=`echo $ref1-$ref2 | bc`; 
                twohigh=`echo $ref1+$ref2*2 | bc` ; twolow=`echo $ref1-$ref2*2 | bc`;
                echo ABOVE AVG; seqlen.awk $file | awk -v HIGH="$high" -v TWOHIGH="$twohigh" '{
                        if ($2>=TWOHIGH) {print $0 "**"}
                        else if ($2>=HIGH) {print $0 "*"}
                        else {print $0}
                        }'
                echo; echo BELOW AVG; seqlen.awk $file | awk -v LOW="$low" -v TWOLOW="$twolow" '{
                        if ($2<=TWOLOW) {print $0 "**"}
                        else if ($2<=LOW) {print $0 "*"}
                        else {print $0}
                        }'         
        done < data.txt ;  
        echo _____________________________ ; 
done

Show how significant the chromosome average repeat counts are using seqlen.awk from my bin and awk programming
for file in *gly*500*tsv; do echo $file; awk 'FNR>1 {print $1}' $file | sort | uniq -c | cut -d ' ' -f1| awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}}
          END {for (i=1;i<=NF;i++) {
          printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
         }' > data.txt ; cat data.txt;
        while IFS=" " read -r ref1 ref2 remainder; do
                high=`echo $ref1+$ref2 | bc` ; low=`echo $ref1-$ref2 | bc`;
                twohigh=`echo "$ref1+$ref2*2" | bc` ; twolow=`echo "$ref1-$ref2*2" | bc`;
                echo ABOVE AVG; awk 'FNR>1 {print $1}' $file | sort | uniq -c | awk -v HIGH="$high" -v TWOHIGH="$twohigh" '{
                        if ($1>=TWOHIGH) {print $0 "**"}
                        else if ($1>=HIGH) {print $0 "*"}
                        else {print $0}
                        }'
                echo; echo BELOW AVG; awk 'FNR>1 {print $1}' $file | sort | uniq -c| awk -v LOW="$low" -v TWOLOW="$twolow" '{
                        if ($1<=TWOLOW) {print $0 "**"}
                        else if ($1<=LOW) {print $0 "*"}
                        else {print $0}
                        }'
        done < data.txt ;
        echo _____________________________ ;
done


Notes about computational biology or related subjects
File formats
"GFF is a standard file format for storing genomic features in a text file. GFF stands for Generic Feature Format. GFF files are plain text, 9 column, tab-delimited files. GFF databases also exist. They use a schema custom built to represent GFF data. GFF is frequently used in GMOD for data exchange and representation of genomic data."
- http://gmod.org/wiki/GFF3

BED format

Python codes
Get currently imported modules
- import sys
- sys.module.keys()

### Clear the console
- with ANSCI commands
  	print('\033c', end='')
- using the terminal command wit os
  	import os
	os.system('cls' if os.name == 'nt' else 'clear')
