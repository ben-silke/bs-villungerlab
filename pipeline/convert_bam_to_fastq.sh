

#!/bin/bash

# Move to the data folder
cd ../data/


# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 24 36 48)
# can use all treatments in the future, but for now we will just used Nutl because of memory issues on the computer.
# treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")
treatments=("Nutl")


# This command should run the samtools fastq command for all treatments/ times.
for treatment in ${treatments[@]}; do
    for time in ${times[@]}; do
        time_folder="${treatment}_${time}.bam"
        folder="/organised/$treatment/$time_folder"
        samtools fastq "${folder}output_merged_${time_folder}.bam" > "${folder}/output_merged_${time_folder}.fastq"
        # see {folder}/
        # samtools fastq "${folder}/output_merged_${time_folder}.bam" > "${folder}/output_merged_${time_folder}.fastq"
    done
done