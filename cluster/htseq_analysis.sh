#!/bin/bash
# Mapping step using STAR

# Move to the data folder
datadir="/nobackup/lab_villunger/bsilke"
cd $datadir

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 20 24 36 48)
replicates=("r1" "r2" "r3" "r4" "r5" "r6")
treatments=("ZM" "Nutl" "Noc" "Etop" "DHCB")

### ATTENTION
treatments=("ZM")

# This command should run the STAR mapping process for all treatments/ times.
for treatment in ${treatments[@]}; do
    output_htseq="$datadir/organised/${treatment}/output_htseq_counts"
    mkdir $output_htseq

    # where do you store stars?
    night_sky="$datadir/organised/${treatment}/output_STAR"
    mkdir $night_sky
    echo "Night Sky: $night_sky"

    input_trimmed_folder="$datadir/organised/${treatment}/output_trimmed"

    for time in ${times[@]}; do
        for replicate in ${replicates[@]}; do
            # /Volumes$datadir/organised/ZM/output_STAR/ZM_24_r4_Aligned.out.sam
            file="${treatment}_${time}_${replicate}_Aligned.out.sam"
            input_file="$datadir/organised/ZM/output_STAR/${file}"
            echo "HTSeq is trying to count in $file"
            output_file="${night_sky}/${file}"

            htseq-count -m union -i gene_id -a 10 --stranded=no $input_file "$datadir/references/Homo_sapiens.GRCh38.109.gtf" > "$output_htseq/htseq_count_${file}.counts"   # output files


        done
    done
done