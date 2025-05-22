#!/bin/bash

# Shell Jutsu: The Way of the Command Line 
# Just some of my command-line jutsu for data wrangling and analysis. Regex throws, awk strikes, and tarball flips—perfect for wrangling massive datasets with minimal motion.
# *"Absorb what is useful, discard what is not, add what is uniquely your own." - Bruce Lee*
# THIS SCRIPT IS NOT EXECTUBLE AS WRITTEN

## File and directory management
#### List files

ls PATH/TO/DIR
 
#### rename a file

mv filename new_filename

#### Add prefix to filenames

for f in * ; do mv -- "$f" "PREFIX_$f" ; done

#### Remove suffix from file names

for f in *; do mv "${f}" "${f/_SUFFIX/}"; done

#### Change directory

cd PATH/TO/DIR 
 
#### Get current directory name only
	
basename "$PWD" 
#OR
basename `pwd` # Get name of current directory only instead of full path

#### Make directory and move into it

mkdir foo && cd "$_" 

#*$_ is a special parameter that holds the last argument of the previous command. Why use double quotes? The quotes around $_ make sure it works even if the folder name contains spaces.*

#### Make directory and copy files into it

mkdir -p /foo/bar && cp myfile.txt $_ 

#### Sync folders and delete or update files

rsync -avu --delete "A/" "B/" # Sync folders and delete or update files

#### Make file executable 

chmod +x file.sh 

#- $#   # Stores the number of command-line arguments that were passed to the shell program.
#- $?   # Stores the exit value of the last command that was executed.
#- $0   # Stores the first word of the entered command (the name of the shell program).
#- $*   # Stores all the arguments that were entered on the command line ($1 $2 ...).
#- "$@"  Stores all the arguments that were entered on the command line, individually quoted ("$1" "$2" ...).

#### Add a folder or file to your PATH variable

vim ~/.bashrc

#### <ins>The following commands are useful for working with compressed files</ins>

#### Create a tar zipped file(s)

tar -czf  name-of-archive.tar.gz /path/to/directory-or-file /path/to/another-directory-or-file …

#### Extract and manipulate zipped files

tar -xvzf filename.tar.gz # Extracts a tarball created by the Tar compression utility

#### Extracts a single, compressed file created by the Gzip compression utility

gunzip filename.gz 

#### View contents of zipped file without extracting

tar -tf filename.tar.gz  

#OR

vim filename.tar.gz # This method allows the user to also see the file contents.

#OR

less filename.tar.gz

#OR

zcat filename.tar.gz

#### To extract lines maching a pattern

zgrep PATTERN  filename.tar.gz 

#### View the difference between files

diff filename1.tar.gz  filename2.tar.gz 

## The *find* command
#### Find files created after date and with a suffix (e.g. '.json')

find . -type f -name "*.json" -newermt "2024-08-01" -exec ls -l {} +

#### Find filed created within a certain date range and matching a pattern

find . -type f -name "*.json" -newermt "2024-09-17" ! -newermt "2024-09-19" -exec ls {} +  

#### Find and remove files matching a PATTERN

find . -maxdepth 1 ! -name "*fna" -exec rm -r {} + # The option maxdepth restricts the command to only run in the current directory

#### Find multiple file patterns at once (could this be simpler?)

find . \( -name "G*glys*gff" -o -name "G*glyd*gff" \)
 
##### I can also *grep* to do this with the syntax 'PATTERN|PATTERN'

#### Get information on running processes

top 

#### Get the location where an executable file is stored. This can be used to check if a program or it's dependencies are installed

which [COMMAND]

## Analyzing and manipulating files

#### Chain linux commands executing subsequent commands only if prior return with zero exit status

Command1 && command2 

#### Run loop on select files

for file in *; do COMMANDS; done

#### Transpose a column to a row

[TABULAR_FEED] |  tr '\n' ' '

#### Transpose tab-delimited line to column

tr '\t' '\n' file 

#### Cut columns except selected

cut -d$'\t' -f 1-10

#### use the word count command to print the number of lines

wc -l 

#### Get lines matching a pattern

grep 'PATTERN' filename 

#### Get blocks of lines matching a pattern using grep

#*Grep has the following options. See the man page for more information:*
#- *-A [num] Print num lines of trailing context after each match. See also the -B and -C options.*
#- *-B [num] Print num lines of leading context before each match. See also the -A and -C options.*
#- *-C [num] Print num lines of leading and trailing context surrounding each match. The default is 2 and is equivalent to -A 2 -B 2. Note, no whitespace may be given between the option and its argument.*

#### Get lines matching one or more patterns and send to a file. Here are severable viable options.

grep ' P \| CA ' file > new_file

grep -E ' (P|CA) ' file > new_file

awk '/ P / || / CA /' file

#- From the man page of grep: -e pattern, --regexp=pattern ; Specify a pattern used during the search of the input: an input line is selected if it matches any of the specified patterns. This option is most useful when multiple -e options are used to specify multiple patterns, or when a pattern begins with a dash (‘-’).

cat file | awk '{print NF}' # Count how many columns in file

awk '!($10="") {print}' filename # Remove a column using awk

awk '$1~/Wm82.gnm2.ann1/ {print}' # Search file for pattern in particular col and print

cat filename.tsv | awk -v max=0 '{if($1>max){want=$2; max=$1}}END{print want}' # Get max numerical value in column
# This works on two column data where you 'want' the maximum value over all other values. 'Max' is a placeholder for the current max value before analyzing the entire dataset.

#### Get the nearest rounded number

var=2.5
echo $var | awk '{print int($1+0.5)}'

#### Use shell variable in awk regex pattern

for x in {01..20}; do 
 	echo Gm$x ; 
   	awk 'FNR>5 && '/'.Gm'$x'/ && $3>=90 {print $9}' out.blastn.Gmr12_rep.glyma.*ISU* | histogram -n | grep -v '0' ; 
     	echo;
done

#### Case-insensitive regex search for negative of pattern at start of pattern

/^[^PATTERN]/I

#### Calculate the average of a column
	
head -15 *Wm82_I*gnm4*.gff | awk '{print $9}'| grep -o -P '(?<=value=).*(?=;ID)'
 
So from the column below: 
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

#I get only the actual values designated by value.
#Now I pipe that to the commands below to get the average.

-> | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }'

#### Calculate the median of a column.

 [TABULAR_FEED]-> | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'
  
#**05-21-25: Check validity of this code; Until then, use an installed program or script go calculate the median :-)**

#**05-21-25:PICK UP EDITING HERE!**
#*Warning: Entries below the following line have not yet been proofread from draft form. Check back soon for more updates on this script!*

----------------------------------------------------------------------------------------------------------------------------------------------------------
#### Calculate the average and standard deviation (population) of a column
awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} 
        END {for (i=1;i<=NF;i++) {
        printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
        }' file.dat >> aver-std.dat

#### Calculate the standard deviation (sample) of a column

awk '{sum+=$0;a[NR]=$0}END{for(i in a)y+=(a[i]-(sum/NR))^2;
 print sqrt(y/(NR-1))}' $filename

# "The two forms of standard deviation are relevant to two different types of variability. One is the variability of values within a set of numbers and one is an estimate of the variability of a population from which a sample of numbers has been drawn.
#The population standard deviation is relevant where the numbers that you have in hand are the entire population, and the sample standard deviation is relevant where the numbers are a sample of a much larger population.
#For any given set of numbers the sample standard deviation is larger than the population standard deviation because there is extra uncertainty involved: the uncertainty that results from sampling. See this for a bit more information: Intuitive explanation for dividing by 𝑛−1 when calculating standard deviation?
#For an example, the population standard deviation of 1,2,3,4,5 is about 1.41 and the sample standard deviation is about 1.58."
# from https://stats.stackexchange.com/questions/485326/confused-when-to-use-population-vs-sample-standard-deviation-in-engineering-test
			    
## Shell arithmetic

echo $((1+1))
	2
 
some_var=2
my_var=$((1+1+some_var))
echo $my_var
  	4

### Help

man cd # Access the manual page for a command

### Execute commands from history
#*see the bash manual page for "history expansion" using 
LESS='+/^HISTORY EXPANSION' man bash*

#In a bash shell...

!636 # run command 636 from history 

!! # run last executed command

!-2 # run the second to last command


cat !$

### Screen Output

echo "" # Print a new line:

echo $'hello\nworld' # print hello world on separate lines
	hello
	world
  
*$'' strings use ANSI C Quoting. The word expands to a string, with backslash-escaped characters replaced as specified by the ANSI C standard.*

printf '%.0s\n' {1..3} # Print multiple blank lines (using bash, ksh93, or zsh at least)

echo {01..20} # list a sequence of numbers

echo this text goes into a file > file.txt

#### Clear the screen

clear

#### Redirect std-error and std-out to a file

[command] 2>&1

# > redirects output to a file, overwriting the file.
# >> redirects output to a file appending the redirected output at the end.
# Standard output is represented in bash with number 1 and standard error is represented with number 2. They are separate, so the user can redirect them to different files.
# 2>&1 redirects the standard error to the standard output so they appear together and can be jointly redirected to a file. (Writing just 2>1 would redirect the standard error to a file called "1", not to standard output.)

### More, less, and most
#See https://askubuntu.com/questions/1191862/what-is-the-difference-between-more-and-less-commands

### FASTA files
#### Split a FASTA file by ID
	
awk '/^>/{if (f) close(f); f=substr($0,2) ".fasta"} {print > f}' input.fasta
#*If sequence IDs have special characters or spaces, consider sanitizing them before using this method.*


# Update a directory from another
rsync -avu --delete "/home/user/A/" "/home/user/B"
# -a Do the sync preserving all filesystem attributes
# -v run verbosely
# -u only copy files with a newer modification time (or size difference if the times are equal)
# --delete delete the files in target folder that do not exist in the source

# Copy multiple files
cp /home/usr/dir/{file1,file2,file3,file4} /home/usr/destination/
#The syntax uses the cp command followed by the path to the directory the desired files are located in with all the files you wish to copy wrapped in brackets and separated by commas.
#Make sure to note that there are no spaces between the files. The last part of the command, /home/usr/destination/, is the directory you wish to copy the files into.
#or if the all the files have the same prefix but different endings you could do something like this:
cp /home/usr/dir/file{1..4} ./ # Where file1,file2,file3 and file4 would be copied.




# Read list items into  loop from text file
for file in `cat test.txt`; do echo $file; done

# Scripting: check existence of dir and variable
[-d "path/to/dir']
[-z ${var+x} ] # if variable does not exist, then null. Otherwise the variable is subbed by x.


# Sort by one column and then the next
sort -k 1,1 -k2,2n file

# Remove duplicates from a fasta file using seqkit
seqkit rmdup -s < in.fa > out.fa

# Remove line breaks from sequences in a fasta file
awk '!/^>/ { printf "%s", $0; n = "\n" } 
  /^>/ { print n $0; n = "" }
  END { printf "%s", n }' input.fasta

  

# Split a string into chunks of 6 and append a prefix to it
echo Hello World | fold -w6 | sed -e 's/^/chunk_/'

# Use desktop calculator command, if installed
dc

# Show how significant the average lengths of sequences in a Fasta file (e.g. set of chromosomes) are using seqlen.awk from my bin and awk programming
for file in gly*/*main.fna; do echo $file; seqlen.awk $file | cut -f2 | awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}}                     
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

# Show how significant the chromosome average repeat counts are using seqlen.awk from my bin and awk programming
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

