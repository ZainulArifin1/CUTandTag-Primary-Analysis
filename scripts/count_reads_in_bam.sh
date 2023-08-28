bam_files=(/home/zarifin_1/cutandtag/output/run64_k27/bam/*.bam)

# Output file
output_file="read_counts_run64_k27.txt"

# Loop through each BAM file and execute the command, appending output to the file
for bam_file in "${bam_files[@]}"; do
    unique_read_count=$(samtools view -F 0x4 "$bam_file" | cut -f 1 | sort | uniq | wc -l)
    echo "Unique read count in $bam_file: $unique_read_count" >> "$output_file"
done