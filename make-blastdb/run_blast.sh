#!/bin/bash

# Path to the input query FASTA file
query_fasta="/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Abundance/FASTA/20250123_nifH_pos_Clostridium_isolates_full16S_con_seq.fa"

# Path to blast database
blast_db="/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Abundance/FASTA/blastdb/blastdb"

# Path to the output directory
output_dir="/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Abundance/FASTA/blast_results"

# Check if query FASTA file exists
if [[ ! -f $query_fasta ]]; then
  echo "Error: Query FASTA file not found!"
  exit 1
fi

# Loop through each sequence in the query FASTA file
awk '/^>/ {if (seq) {print header"\n"seq}; header=$0; seq=""} !/^>/ {seq=seq$0} END {if (seq) {print header"\n"seq}}' "$query_fasta" | \
while read -r header; do
  # Read the corresponding sequence
  read -r seq

  # Extract the sequence name from the header (remove the ">" symbol)
  seq_name=$(echo "$header" | sed 's/^>//' | tr -d ' \t')

  # Create a unique temporary FASTA file for the current sequence
  temp_query=$(mktemp /tmp/temp_XXXXXX.fasta)

  # Write the current sequence to the temporary FASTA file
  echo -e "$header\n$seq" > "$temp_query"

  # Run the BLAST search and save the output to a file named after the sequence
  blastn -query "$temp_query" -db "$blast_db" -outfmt "6 qseqid sseqid pident length evalue bitscore" -evalue 1e-5 | \
  tee >(awk 'BEGIN {print "qseqid\tsseqid\tpident\tlength\tevalue\tbitscore"} {print}' > "$output_dir/${seq_name}_blast_results.txt") | \
  awk 'BEGIN {print "qseqid\tsseqid\tpident\tlength\tevalue\tbitscore"} $3 >= 97 {print}' > "$output_dir/filtered_results/${seq_name}_filtered_results.txt"

  # Clean up the temporary file
  rm -f "$temp_query"

  echo "Processed sequence: $seq_name"
done
