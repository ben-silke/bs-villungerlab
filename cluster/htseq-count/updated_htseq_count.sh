#!/bin/bash
# Counting htseq-counts

treatment=$1
echo "We are now operating on $treatment"

# Move to the data folder
datadir="/nobackup/lab_villunger/bsilke"
cd $datadir

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 20 24 36 48)
replicates=("r1" "r2" "r3" "r4" "r5" "r6")
treatments=("ZM" "Nutl" "Noc" "Etop" "DHCB")

output_htseq_count="/nobackup/lab_villunger/bsilke/organised/${treatment}/output_htseq_counts"
mkdir $output_htseq_count

# where do you store stars?
night_sky="/nobackup/lab_villunger/bsilke/organised/${treatment}/output_STAR"
files=()

for time in ${times[@]}; do
    for replicate in ${replicates[@]}; do
        file="${treatment}_${time}_${replicate}"
        echo "Mapping STARS in $file"
        input_file="${night_sky}/${file}Aligned.out.sam"
        files+=("$input_file")
    done
done

echo "All Files: ${files[@]}"
# [bsilke@d002 output_STAR]$ htseq-count -i "gene_id" ZM_0_r1Aligned.out.sam /nobackup/lab_villunger/bsilke/references/star/Homo_sapiens.GRCh38.110.gtf
htseq-count \
-i "gene_id" \
# -n 8 \
# --with-header \
-c "${output_htseq_count}/${treatment}_htseqcount" \
${files[@]} \
"/nobackup/lab_villunger/bsilke/references/star/Homo_sapiens.GRCh38.110.gtf"

# done