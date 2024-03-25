#!/bin/bash

echo "Downloading libraries"	
yes | sudo apt install fastqc
yes | sudo apt install rna-star
yes | sudo apt install samtools
yes | sudo apt install subread
yes | sudo apt install trimmomatic
yes | sudo apt install STAR
yes | sudo apt install bwa
yes | sudo apt install varscan
yes | sudo apt install bedtools
yes | conda install trimmomatic
#yes | conda install STAR
#yes | conda install samtools
#yes | conda install subread
#yes | conda install fastqc
#yes | conda install rna-star
#yes | conda install bwa
#yes | conda install varscan
#yes | conda install bedtools

#find / -type f -name "trimmomatic" 2>/dev/null
