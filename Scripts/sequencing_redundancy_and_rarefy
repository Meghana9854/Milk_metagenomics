mkdir -p ./5.nonpareil/2.redundancy_kmer

##### ESTIMATING SEQUENCING REDUNDANCY WITH NONPAREIL #####
module load nonpareil/3.4.1
# estimating the sequencing effort with the estimated diversity to look at the redundancy
for i in ./5.nonpareil/1.fastq_reads/new_files/*_1.fastq
do
 nonpareil -s "$i" -T kmer -k 20 -X 10000 -t 16 -f fastq -b ./5.nonpareil/2.redundancy_kmer/"$(basename "$i" .fastq)"
done

# moving the npo files needed for downstream analysis into a single folder
mkdir -p ./5.nonpareil/2.redundancy_kmer/npo_files
mv ./5.nonpareil/2.redundancy_kmer/*.npo ./5.nonpareil/2.redundancy_kmer/npo_files

##### RAREFY READS #####
module load seqtk/1.3

mkdir -p ./7.rarefy

# Output directory for subsampled reads
output_dir="./7.rarefy/1.fastq"
mkdir -p "$output_dir"

# change seed for every subsamples (5, 11, ) so you get different subsamples #
# Seeds for subsampling
seeds=(5 11 20 45 62 75 90 135 150 189)

# Iterate over each filename in filenames.txt
while IFS= read -r i; do
    # Reset the subsampling counter for each new file
    sub=0

    # Iterate over each seed
    for seed in "${seeds[@]}"; do
        ((sub++)) # Increment subsampling counter
        
        echo "$i - executing seqtk for seed: $seed and sub: $sub"
        
        # Paths for the original R1 and R2 files
        r1="./3.host_removal/${i}_microbial_R1.fastq.gz"
        r2="./3.host_removal/${i}_microbial_R2.fastq.gz"
        
        # Output paths for the subsampled R1 and corresponding R2 files
        out_r1="${output_dir}/${i}_sub${sub}_R1.fastq"
        out_r2="${output_dir}/${i}_sub${sub}_R2.fastq"
        
        # make a temporary fastq file that contains the subsample of R1 
        temp_r1="${output_dir}/temp_${i}_sub${sub}_R1.fastq"
        
        # Subsample R1
        seqtk sample -s"$seed" "$r1" 1000000 > "$temp_r1"
        
        # Extract read IDs from the subsampled R1, removing the /1 suffix
        awk 'NR % 4 == 1' "$temp_r1" > "${output_dir}/${i}_sub${sub}_ids_for_r1.txt"
        
        # Use zgrep to filter R1 reads based on the IDs taken from the temp subsampled R1 file
        zgrep -A 3 -F -f "${output_dir}/${i}_sub${sub}_ids_for_r1.txt" "$r1" | grep -v "^--$" > "$out_r1"
        
        # For R2, append /2 to each ID before filtering
        awk '{sub(/\/1$/, "/2"); print}' "${output_dir}/${i}_sub${sub}_ids_for_r1.txt" > "${output_dir}/${i}_sub${sub}_ids_for_r2.txt"
        
        # Use zgrep to filter R2 reads based on the adjusted IDs
        zgrep -A 3 -F -f "${output_dir}/${i}_sub${sub}_ids_for_r2.txt" "$r2" | grep -v "^--$" > "$out_r2"

    done
done < filenames.txt

module unload seqtk/1.3
