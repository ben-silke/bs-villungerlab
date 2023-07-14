#!/bin/bash
# Mapping step using STAR

# Move to the data folder
datadir="/nobackup/lab_villunger/bsilke"
cd $external_drive

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 20 24 36 48)
replicates=("r1" "r2" "r3" "r4" "r5" "r6")
treatments=("ZM" "Nutl" "Noc" "Etop" "DHCB")

# STAR \
# --runThreadN 2 \
# --genomeDir ../../data/references/STAR \
# --readFilesIn "specific to the files"

ls
# This command should run the STAR mapping process for all treatments/ times.
for treatment in ${treatments[@]}; do
    output_htseq="$datadir/organised/${treatment}/output_htseq_counts"
    mkdir $output_htseq

    # where do you store stars?
    night_sky="$datadir/organised/${treatment}/output_STAR"
    mkdir $night_sky
    

    input_trimmed_folder="$datadir/organised/${treatment}/output_trimmed"

    for time in ${times[@]}; do
        for replicate in ${replicates[@]}; do
            file="${treatment}_${time}_${replicate}"

            echo "Mapping STARS in $file"

            input_trimmed_file="${input_trimmed_folder}/trimmed_fq_${file}.fastq"
            output_file="${night_sky}/${file}"

            STAR \
            --runThreadN 4 \
            --genomeDir $datadir/star_index \
            --readFilesIn "$input_trimmed_file" \
            --outFileNamePrefix "${night_sky}/${file}_"
            # Now we should delete the old file
            # rm $input_trimmed_file
        done
    done
done