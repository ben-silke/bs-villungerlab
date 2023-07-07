#!/bin/bash

# Move to the data folder
# cd ../data/
main_path=$1
cd $main_path

# Move to the data folder
# external_drive="/Volumes/bs_external/villunger"
# cd $external_drive

ls
# Include all times, can exclude based upon experiement later.
# TODO: make these env variables so that you only need to update this once
times=(0 8 12 16 20 24 36 48)
# can use all treatments in the future, but for now we will just used Nutl because of memory issues on the computer.
treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")
replicates=("r1" "r2" "r3" "r4" "r5" "r6")

# SPECIFIC TREATMENT
# treatments=("ZM")
# ZM does not have 8 or 12 so can exclude
# times=(0 16 24 36 48)
treatments=("ZM" "Nutl" "Noc" "Etop" "DHCB")

# This command should run the samtools fastq command for all treatments/ times.
for treatment in ${treatments[@]}; do
    output_fastq_folder="organised/${treatment}/output_merged_fq"
    mkdir $output_fastq_folder
    for time in ${times[@]}; do
        for replicate in ${replicates[@]}; do

            file="${treatment}_${time}_${replicate}"
            echo "Converting $file"
            output_merged_folder="organised/${treatment}/output_merged"

            samtools fastq "$output_merged_folder/output_merged_${file}.bam" > "$output_fastq_folder/output_merged_fq_${file}.fastq"

            # We should probably remove the bam file now to make space
            # rm "organised/${treatment}/output_merged_${file}.bam"
        done
    done
done