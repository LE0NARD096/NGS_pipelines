#!/bin/bash

data_dir="DataTP1"
results_dir="ResultsTP1"

if [ ! -d "$data_dir" ]; then
    mkdir "$data_dir"
    echo "Downloading files and setting up directory..."

    # Download the files (force download if they already exist)
    wget -N -P "$data_dir" http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz
    wget -N -P "$data_dir" http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz
    wget -N -P "$data_dir" ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz

    echo "Extracting files"
    tar -zxvf "$data_dir/TPrnaseq.tar.gz" -C "$data_dir"
    gunzip "$data_dir/chr18.fa.gz"
    gunzip "$data_dir/gencode.v24lift37.basic.annotation.gtf.gz"

    echo "Directory and files are set up."
else
    echo "The '$data_dir' directory already exists. Skipping download."
fi

# STAR manual

if [ ! -d GenDir1 ]; then
        echo "Downloading files and setting up directory..."
        mkdir GenDir1

         STAR --runMode genomeGenerate --runThreadN 4 \
        --genomeDir GenDir1 \
        --genomeFastaFiles $data_dir/chr18.fa \
        --sjdbGTFfile $data_dir/gencode.v24lift37.basic.annotation.gtf

  echo "Directory and files are set up."
else
  echo "The 'GenDir1' directory already exists. Skipping download."
fi

# Trimming and adapter removal

echo "Trimming files..."

if [ ! -d "$results_dir" ]; then
    mkdir "$results_dir"
fi

# Utilizza un ciclo for per iterare su tutti i file .fastq nella directory Data
for R1_input in "$data_dir"/*.fastq; do
    if [ -e "$R1_input" ]; then
        R2_input="${R1_input/R1/R2}"  # Trova il corrispondente file R2
        output_base=$(basename "${R1_input}" .fastq)
        /var/lib/miniforge/share/trimmomatic-0.39-2/trimmomatic PE "$R1_input" "$R2_input" -baseout "$results_dir/${output_base}.fastq" LEADING:20 TRAILING:20 MINLEN:50
    fi
done

# Mapping

echo "Reading Genome..."
if [ ! -d "$results_dir/STAR_FILES" ]; then
    mkdir "$results_dir/STAR_FILES"
fi

for file in "$data_dir"/*R1.fastq; do
    r1_file="$file"
    r2_file="${r1_file/R1.fastq/R2.fastq}"

    echo $r1_file

    if [ -e "$r1_file" ] && [ -e "$r2_file" ]; then
        sample_name=$(basename "$r1_file" | cut -d. -f1)

        STAR --runThreadN 4 --outFilterMultimapNmax 1 \
            --genomeDir GenDir1 \
            --outSAMattributes All --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix "$results_dir/STAR_FILES/$sample_name" \
            --readFilesIn "$r1_file" "$r2_file"

        samtools index "$results_dir/STAR_FILES/${sample_name}Aligned.sortedByCoord.out.bam"
    fi
done

echo "Counting..."

featureCounts -p -t exon -g gene_id -a "$data_dir/gencode.v24lift37.basic.annotation.gtf" -o "$results_dir/count.out" "$results_dir/STAR_FILES"/*Aligned.sortedByCoord.out.bam

perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' "$data_dir/gencode.v24lift37.basic.annotation.gtf" | sort | uniq > "$results_dir/encode-to-hugo.tab"
sort "$results_dir/count.out" > "$results_dir/temp1"
sort "$results_dir/encode-to-hugo.tab" > "$results_dir/temp2"
join "$results_dir/temp1" "$results_dir/temp2" | grep "chr18" > "$results_dir/temp3"
awk '{print $NF,$(NF-7), $(NF-6), $(NF-5), $(NF-4), $(NF-3), $(NF-2), $(NF-1)}' "$results_dir/temp3" > "$results_dir/count.txt"

echo "Script completed"
