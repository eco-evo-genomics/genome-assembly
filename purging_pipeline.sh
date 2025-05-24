#!/bin/bash -l

# =============================================================================
# Purge_Dups Genome Refinement Pipeline: Primary and Haplotig Cleanup
#
# This script refines a genome assembly using PacBio reads and `purge_dups`,
# removing haplotypic duplications from both the primary assembly and a
# merged haplotig assembly.
#
# Requirements:
# - `purge_dups` compiled binaries in your PATH (or adjust paths below)
# - `minimap2` installed and accessible
#
# Notes:
# - Outputs are written safely with unique names for primary and haplotig runs
# - Script uses multithreading (36 CPUs)
#
# =============================================================================

# -----------------------------
# User-defined variables
# -----------------------------

# Primary assembly (from hifiasm or other assembler)
PRI_ASM="primary_assembly.fa"

# Alternative assembly (usually the a_ctg.fa file from hifiasm)
ALT_ASM="alternate_assembly.fa"

# PacBio reads (FASTA or FASTQ)
PB_READS="pacbio_reads.fasta"

# Number of threads to use
THREADS=36

# Path to purge_dups binaries
BIN_DIR="bin"

# -----------------------------
# Safety settings
# -----------------------------
set -euo pipefail

# -----------------------------
# Step 1: Align PacBio reads to primary assembly
# -----------------------------
echo "Aligning PacBio reads to primary assembly..."
minimap2 -t $THREADS -xasm20 $PRI_ASM $PB_READS | gzip -c - > pb_reads.paf.gz

# Check if bin directory exists
[ -d $BIN_DIR ] || mkdir -p $BIN_DIR

echo "Generating coverage stats..."
$BIN_DIR/pbcstat pb_reads.paf.gz

echo "Calculating coverage cutoffs..."
$BIN_DIR/calcuts PB.stat > cutoffs 2> calcuts.log

echo "Splitting primary assembly..."
$BIN_DIR/split_fa $PRI_ASM > ${PRI_ASM}.split

echo "Self-aligning primary assembly..."
minimap2 -t $THREADS -xasm5 -DP ${PRI_ASM}.split ${PRI_ASM}.split | gzip -c - > ${PRI_ASM}.split.self.paf.gz

echo "Purging duplicates from primary assembly..."
$BIN_DIR/purge_dups -2 -T cutoffs -c PB.base.cov ${PRI_ASM}.split.self.paf.gz > dups.bed 2> purge_dups.log

echo "Extracting cleaned primary and haplotig sequences..."
$BIN_DIR/get_seqs -e dups.bed $PRI_ASM

# At this point:
# - `purged.fa` contains the cleaned primary assembly
# - `hap.fa` contains putative haplotigs

# -----------------------------
# Step 2: Merge and refine haplotigs
# -----------------------------
echo "Merging hap.fa with $ALT_ASM to create merged_hap.fa..."
cat hap.fa $ALT_ASM > merged_hap.fa

echo "Aligning PacBio reads to merged haplotigs..."
minimap2 -t $THREADS -xasm20 merged_hap.fa $PB_READS | gzip -c - > pb_reads.merged.paf.gz

echo "Generating coverage stats for merged haplotigs..."
$BIN_DIR/pbcstat pb_reads.merged.paf.gz

echo "Calculating cutoffs for merged haplotigs..."
$BIN_DIR/calcuts PB.stat > cutoffs.merged 2> calcuts_merged.log

echo "Splitting merged haplotigs..."
$BIN_DIR/split_fa merged_hap.fa > merged_hap.fa.split

echo "Self-aligning merged haplotigs..."
minimap2 -t $THREADS -xasm5 -DP merged_hap.fa.split merged_hap.fa.split | gzip -c - > merged_hap.fa.split.self.paf.gz

echo "Purging duplicates from merged haplotigs..."
$BIN_DIR/purge_dups -2 -T cutoffs.merged -c PB.base.cov merged_hap.fa.split.self.paf.gz > dups_merged.bed 2> purge_dups_merged.log

echo "Extracting final cleaned haplotig set..."
$BIN_DIR/get_seqs -e dups_merged.bed merged_hap.fa

# -----------------------------
# Done
# -----------------------------
echo "Purge_Dups pipeline complete."
echo "Final outputs:"
echo " - purged.fa (cleaned haplotig assembly)"
echo " - hap.fa (leftover redundant sequences)"
