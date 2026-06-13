#!/bin/bash

#SBATCH --job-name=con_gen_05
#SBATCH --output=con_gen_05_%j.out
#SBATCH --error=con_gen_05_%j.err
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=48G

# Load necessary modules
module load tools  # Ensure this loads bedtools and other required tools

# Define variables
REFERENCE="/path/to/reference/genome.fasta"
INPUT_BAM="/path/to/alignment/alignment.sorted.bam"
MIN_COVERAGE=100
MIN_ALLELE_FREQ=0.5  # Minimum allele frequency threshold (e.g., 50% = 0.5)
MIN_QUAL=20  # Minimum base quality threshold for bcftools mpileup
OUTPUT_DIR="/path/to/output/consensus"
GENOME_FILE="$OUTPUT_DIR/genome.txt"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR || { echo "Error: Failed to create output directory $OUTPUT_DIR"; exit 1; }

# Generate genome file for bedtools
cut -f1,2 ${REFERENCE}.fai > $GENOME_FILE

# Generate the regions.bed file based on minimum coverage
bedtools genomecov -ibam "$INPUT_BAM" -bga |
    awk -v min_cov=$MIN_COVERAGE '$4 >= min_cov {print $1"\t"$2"\t"$3}' |
    bedtools merge -i - > "$OUTPUT_DIR/regions.bed" || { echo "Error running bedtools genomecov."; exit 1; }

# Check if regions.bed file is not empty
if [ ! -s "$OUTPUT_DIR/regions.bed" ]; then
    echo "Error: regions.bed file is empty or not created."
    exit 1
fi

# Display the contents of regions.bed
echo "Contents of regions.bed:"
cat "$OUTPUT_DIR/regions.bed"

# Generate the VCF file with variants (Ensure mixed genotypes are retained)
bcftools mpileup -Ou -Q $MIN_QUAL -f "$REFERENCE" -a FORMAT/DP,AD "$INPUT_BAM" |
    bcftools call -m -Oz -o "$OUTPUT_DIR/variants.vcf.gz" || { echo "Error generating VCF file."; exit 1; }

# Filter the VCF file to **retain** heterozygous sites for IUPAC encoding
bcftools filter -e "FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1]) < $MIN_ALLELE_FREQ" -Oz \
-o "$OUTPUT_DIR/filtered_variants.vcf.gz" "$OUTPUT_DIR/variants.vcf.gz" \
 || { echo "Error filtering VCF file."; exit 1; }

# Normalize VCF to merge alleles and ensure correct IUPAC encoding
bcftools norm -m -both -f "$REFERENCE" "$OUTPUT_DIR/filtered_variants.vcf.gz" -Oz -o "$OUTPUT_DIR/normalized_variants.vcf.gz"

# Index the normalized VCF file
bcftools index "$OUTPUT_DIR/normalized_variants.vcf.gz"

# Generate consensus sequence using IUPAC ambiguity codes
bcftools consensus -f "$REFERENCE" -o "$OUTPUT_DIR/consensus.fa" -I "$OUTPUT_DIR/normalized_variants.vcf.gz" || { echo "Error generating consensus sequence."; exit 1; }

echo "Consensus sequence generated at $OUTPUT_DIR/consensus.fa"
