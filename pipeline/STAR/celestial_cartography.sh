#!/bin/bash
# Mapping step using STAR

# Move to the data folder
cd ../../data/

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 24 36 48)
replicates=("r1" "r2" "r3" "r4" "r5" "r6")
treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")

# SPECIFIC TREATMENT
treatments=("ZM")
# ZM does not have 8 or 12 so can exclude
times=(0 16 24 36 48)

# STAR \
# --runThreadN 2 \
# --genomeDir ../../data/references/STAR \
# --readFilesIn "specific to the files"


# This command should run the STAR mapping process for all treatments/ times.
for treatment in ${treatments[@]}; do
    for time in ${times[@]}; do
        for replicate in ${replicates[@]}; do

            file="${treatment}_${time}_${replicate}"

            echo "Mapping STARS in $file"

            # folder for salmon outputs == aquarium obviously
            night_sky="organised/${treatment}/output_salmon"
            mkdir $night_sky

            input_trimmed_folder="organised/${treatment}/output_trimmed"
            input_trimmed_file="${input_trimmed_folder}/trimmed_fq_${file}.fastq"

            STAR \
            --runThreadN 2 \
            --genomeDir references/STAR \
            --readFilesIn "$input_trimmed_file" \
            --outFileNamePrefix "$night_sky"

            # Now we should delete the old file
            # rm $input_trimmed_file
        done
    done
done