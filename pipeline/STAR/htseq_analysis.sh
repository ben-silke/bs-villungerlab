

#!/bin/bash

# Move to the data folder
external_drive="/Volumes/bs_external/villunger"
cd $external_drive

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 20 24 36 48)
replicates=("r1" "r2" "r3" "r4" "r5" "r6")
treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")


# This command should run the STAR mapping process for all treatments/ times.
for treatment in ${treatments[@]}; do
    output_htseq="organised/${treatment}/output_htseq_counts"
    mkdir $output_htseq

    # where do you store stars?
    night_sky="organised/${treatment}/ouput_STAR"
    mkdir $night_sky

    input_trimmed_folder="organised/${treatment}/output_trimmed"

    for time in ${times[@]}; do
        for replicate in ${replicates[@]}; do

            file="${treatment}_${time}_${replicate}"

            # if htseq is not installed on the cluster

            htseq-count \
            -m union \
            -i gene_name \
            -a 10 \
            --stranded=no \
            $night_sky \              # aligned out files
            "$external_drive/STAR/ensemble_annotation/homo_sapiens/Homo_sapiens.GRCh38.109.gff3" \  # annotation file
            "$output_htseq/htseq_count_${file}.counts"   # output files
        done
    done
done