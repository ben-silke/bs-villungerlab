#!/bin/bash

# Move to the data folder
cd ../data/

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 24 36 48)
replicates=("r1" "r2" "r3" "r4" "r5" "r6")
treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")

# SPECIFIC TREATMENT
treatments=("ZM")
# ZM does not have 8 or 12 so can exclude
times=(0 16 24 36 48)

# salmon quant 
# -i salmon_index/Hsapiens_index_k17 
# -l A 
# -r trimmed/RPE_Nutl_t0_r1_tc.fastq 
# -p 8 
# --validateMappings 
# -o quantification/RPE_Nutl_t0_r1   

# Build index if you have not already built the index.
# salmon index -t references/human_index_salmon/gentrome.fa -i references/human_index_salmon/Hsapiens_index_k17 --decoys references/human_index_salmon/decoys.txt -k 17 

# This command should run the salmon quant command for all treatments/ times.
for treatment in ${treatments[@]}; do
    for time in ${times[@]}; do
        for replicate in ${replicates[@]}; do

            file="${treatment}_${time}_${replicate}"
            folder="organised/$treatment"

            echo "Fishing in $file"

            # folder for salmon outputs == aquarium obviously
            aquarium="organised/${treatment}/output_salmon"
            mkdir $aquarium

            input_trimmed_folder="organised/${treatment}/output_trimmed"
            input_trimmed_file="${input_trimmed_folder}/trimmed_fq_${file}.fastq"

            salmon quant -i "references/human_index_salmon/Hsapiens_index_k17"  -o "$aquarium/salmon_quant_${file}" \
            -l A \
            -r $input_trimmed_file \
            -p 8

            # Now we should delete the old file
            rm $input_trimmed_file
        done
    done
done