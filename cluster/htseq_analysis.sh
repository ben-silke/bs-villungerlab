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

ls
# This command should run the STAR mapping process for all treatments/ times.
for treatment in ${treatments[@]}; do
    output_htseq="/nobackup/lab_villunger/bsilke/organised/${treatment}/output_htseq_counts"
    mkdir $output_htseq

    # where do you store stars?
    night_sky="/nobackup/lab_villunger/bsilke/organised/${treatment}/output_STAR"
    mkdir $night_sky
    echo "Night Sky: $night_sky"

    input_trimmed_folder="/nobackup/lab_villunger/bsilke/organised/${treatment}/output_trimmed"

    for time in ${times[@]}; do
        for replicate in ${replicates[@]}; do
            # /Volumes/nobackup/lab_villunger/bsilke/organised/ZM/output_STAR/ZM_24_r4_Aligned.out.sam
            file="${treatment}_${time}_${replicate}_Aligned.out.sam"

            echo "Mapping STARS in $file"
            output_file="${night_sky}/${file}"

            # --outFileNamePrefix "${night_sky}/${file}"

            # htseq-count [options] <alignment_files> <gtf_file>
            htseq-count \
            $file

            htseq-count \
            -m union \
            -i gene_name \
            -a 10 \
            --stranded=no \
            $night_sky \              # aligned out files
            "/nobackup/lab_villunger/bsilke/references/STAR/ensemble_annotation/homo_sapiens/Homo_sapiens.GRCh38.109.gff3" \  # annotation file
            "$output_htseq/htseq_count_${file}.counts"   # output files


        done
    done
done