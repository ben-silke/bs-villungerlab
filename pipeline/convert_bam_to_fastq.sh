#!/bin/bash

# Move to the data folder
cd ../data/

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 24 36 48)
# can use all treatments in the future, but for now we will just used Nutl because of memory issues on the computer.
# treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")
treatments=("Nutl")
replicates=("r1","r2","r3")

# This command should run the samtools fastq command for all treatments/ times.
for treatment in ${treatments[@]}; do
    for time in ${times[@]}; do
        time_folder="${treatment}_${time}"
        echo "Converting $time_folder"
        output_merged_folder="/organised/${treatment}/output_merged"
        output_fastq_folder="/organised/${treatment}/output_merged_fq"
        mkdir $output_fastq_folder
        samtools fastq "$output_merged_folder/output_merged_${time_folder}.bam" > "$output_fastq_folder/output_merged_fq_${time_folder}.fastq"

        # We should probably remove the bam file now to make space
        # rm "organised/${treatment}/output_merged_${time_folder}.bam"
    done
done