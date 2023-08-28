# Calculate the total size of hg38
total_genome_size=2913022398

# Calculate the total size of the blacklisted regions
blacklisted_size=$(awk '{sum+=$3-$2} END {print sum}' ../data/blacklist/hg38-blacklist.v2.bed)

# Calculate the effective genome size
effective_genome_size=$((total_genome_size - blacklisted_size))

# Print the effective genome size
echo "Effective Genome Size for hg38: $effective_genome_size"
