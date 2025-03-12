run_blast.sh is a script to extract the sequences from a fasta file and run them against a created blast database. 
It then outputs both a filtered and unfiltered directory of the results in the specified output path.
You can change the stringency of the percent identity filter by editing this line of code:
awk 'BEGIN {print "qseqid\tsseqid\tpident\tlength\tevalue\tbitscore"} $3 >= 97 {print}' > "$output_dir/filtered_results/${seq_name}_filtered_results.txt"


