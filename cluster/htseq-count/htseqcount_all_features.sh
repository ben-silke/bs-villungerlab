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

#!/bin/bash

if [[ "$treatment" == "Nutl" ]] || [[ "$treatment" == "Etop" ]]
then
  times=(0 8 12 16 24 48)
else
  times=(0 16 20 24 36 48)
fi

# Print the times array
echo "${times[@]}"

output_htseq_counts="/nobackup/lab_villunger/bsilke/organised/${treatment}/output_htseq_counts_2"
mkdir $output_htseq_counts

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
output_file="${output_htseq_counts}/all_${treatment}_fc.tsv"

htseq-count \
--nonunique all \
--stranded no \
-i "gene_id" \
-n 8 \
--with-header \
-c $output_file \
${files[@]} \
"/nobackup/lab_villunger/bsilke/references/star/Homo_sapiens.GRCh38.110.gtf" 
