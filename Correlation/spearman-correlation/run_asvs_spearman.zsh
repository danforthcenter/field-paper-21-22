# Input arguments
input_file="$1"         # BLAST results file
input_asv_table="$2"    # ASV abundance table (CSV)
output_directory="$3"   # Directory to save results
rscript_path="/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/21_22_R_project/Correlation/spearman_rscript.R"

# Check if input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file not found: $input_file"
    exit 1
fi

# Check if ASV table exists
if [[ ! -f "$input_asv_table" ]]; then
    echo "Error: ASV abundance table not found: $input_asv_table"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_directory"

# Extract unique ASVs from the second column ($sseqid) of the BLAST results
asvs=($(awk -F'\t' 'NR>1 {print $2}' "$input_file" | sort -u))
echo "Extracted ASVs:"
echo "$asvs"

# Check if any ASVs were found
if [[ ${#asvs[@]} -eq 0 ]]; then
    echo "Error: No ASVs found in input file."
    exit 1
fi

# Loop through each ASV and call the R script
echo "Processing ASVs..."
for asv in "${asvs[@]}"; do
    if [[ -n "$asv" ]]; then
        echo "🔹 Processing ASV: $asv"
        echo "Calling R script with arguments: ASV=$asv, Table=$input_asv_table, Output=$output_directory"
        echo "Rscript Path: $rscript_path"

        # Call the R script
        Rscript "$rscript_path" "$asv" "$input_asv_table" "$output_directory"
	exit_code=$?
        # Check the exit status of the R script
        if [[ $exit_code -ne 0 ]]; then
            echo "Error: R script failed for ASV: $asv. Stopping further execution."
            exit 1  # Stop the script if the R script fails
        fi

        echo "✔ Completed ASV: $asv"
    fi
done

