#!/bin/bash

# Move to the data folder
cd ../data/


# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 24 36 48)
# can use all treatments in the future, but for now we will just used Nutl because of memory issues on the computer.
# treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")
treatments=("Nutl")


# This command should run the samtools merge command for all treatments/ times.
for treatment in ${treatments[@]}; do
    for time in ${times[@]}; do
        time_folder="${treatment}_${time}"
        folder="organised/$treatment/$time_folder"
        echo "Merging $folder"
        # NOTE: We may need to look at sorting these files.
        samtools merge -o "${folder}/output_merged_${time_folder}.bam" $folder/*.bam

        # To remove the files
        # These files are copied from data/raw so can easily be recopied.
        rm $folder/*.bam
    done
done