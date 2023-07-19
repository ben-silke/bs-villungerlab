#!/bin/bash
# Counting Features

treatment=$1
echo "We are now operating on $treatment"

# Move to the data folder
datadir="/nobackup/lab_villunger/bsilke"
cd $datadir

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 20 24 36 48)
replicates=("r1" "r2" "r3" "r4" "r5" "r6")
treatments=("ZM" "Nutl" "Noc" "Etop" "DHCB")

ls
output_feature_counts="/nobackup/lab_villunger/bsilke/organised/${treatment}/output_feature_counts_gencode"
mkdir $output_feature_counts

# where do you store stars?
night_sky="/nobackup/lab_villunger/bsilke/organised/${treatment}/output_STAR_gencode"
files=()

for time in ${times[@]}; do
    for replicate in ${replicates[@]}; do
        file="${treatment}_${time}_${replicate}"

        echo "Mapping STARS in $file"

        input_file="${night_sky}/gc_${file}Aligned.out.sam"
        files+=("$input_file")
    done
done

echo "All Files: ${files[@]}"

output_file="${output_feature_counts}/allgc_${file}_fc"

featureCounts -F GTF \
-a /nobackup/lab_villunger/bsilke/references/star/gencode.v44.annotation.gtf \
-o $output_file \
${files[@]} \
# done