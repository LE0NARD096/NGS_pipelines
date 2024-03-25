#!/bin/bash

data_dir="DataTP2"
result_dir="ResultsTP2"

# Funzione per gestire gli errori
handle_error() {
    echo "Error: $1"
    exit 1
}

# Creazione delle directory se non esistono
mkdir -p "$data_dir" || handle_error "Impossible create directory $data_dir"
mkdir -p "$result_dir" || handle_error "Impossible create directory $result_dir"

echo "Downloading files and setting up directory..."
# Scarica i file usando i cookie temporanei
file_id="1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe"
cookie=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate "https://docs.google.com/uc?export=download&id=$file_id" -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$cookie&id=$file_id" -O "$data_dir/patient7.tar.gz" && rm -rf /tmp/cookies.txt

# Scarica altri file necessari
wget -N -P "$data_dir" "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz"
wget -N -P "$data_dir" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.annotation.gtf.gz"

# Estrai e decomprimi i file
tar -xzf "$data_dir/patient7.tar.gz" -C "$data_dir" || handle_error "Impossible extract  file patient7.tar.gz"
gunzip "$data_dir/patient7.exome"/* || handle_error "Error during decompression patient7.exome"
gunzip "$data_dir/chr16.fa.gz" || handle_error "Error during decompression chr16.fa.gz"
gunzip "$data_dir/gencode.v24lift37.annotation.gtf.gz" || handle_error "Error during decompression gencode.v24lift37.annotation.gtf.gz"

echo "Indexing chr16..."
bwa index -a bwtsw "$data_dir/chr16.fa" || handle_error "Error during indexing chr16.fa"

echo "Trimming and applying BWA to Norm and Canc genome..."
for F in N T; do
    # Trimming
    /var/lib/miniforge/share/trimmomatic-0.39-2/trimmomatic PE "$data_dir/patient7.exome/TCRBOA7-${F}-WEX-chr16_r1F.fastq" "$data_dir/patient7.exome/TCRBOA7-${F}-WEX-chr16_r2F.fastq" -baseout "$data_dir/${F}.fastq" LEADING:20 TRAILING:20 MINLEN:50

    # BWA
    bwa mem -M -t 2 -A 2 -E 1 "$data_dir/chr16.fa" "$data_dir/${F}_1P.fastq" "$data_dir/${F}_2P.fastq" > "$result_dir/output_${F}.sam"

    # Mapping
    samtools view -S -b "$result_dir/output_${F}.sam" > "$result_dir/output_${F}.bam"
	samtools flagstat "$result_dir/output_${F}.bam" > "$result_dir/flagstat_output_${F}.txt"
    samtools sort "$result_dir/output_${F}.bam" > "$result_dir/sorted_output_${F}.bam"
	samtools index "$result_dir/sorted_output_${F}.bam"
    samtools mpileup -B -A -f "$data_dir/chr16.fa" "$result_dir/sorted_output_${F}.bam" > "$result_dir/${F}_pileup"
done

echo "Calling somatic variants..."
varscan somatic "$result_dir/N_pileup" "$result_dir/T_pileup" "$result_dir/variant_output"

echo "Extracting all somatic mutations from VCF files (SNP and INDEL) and convert to BED format..."
for z in indel snp; do
    grep -i 'somatic' "$result_dir/variant_output.${z}" > "$result_dir/filtered_${z}"
    awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' "$result_dir/filtered_${z}" > "$result_dir/bed_${z}"
    bedtools intersect -a "$data_dir/gencode.v24lift37.annotation.gtf" -b "$result_dir/bed_${z}" > "$result_dir/intersect_${z}"
    grep '\sgene\s' "$result_dir/intersect_${z}" | awk '{print " " $1 " " $4 " " $5 " " $16}' > "$result_dir/results_${z}"
done

echo "Script completed successfully, consult resusts about mutations on ResultsTP2/results_indel and results_snp"
