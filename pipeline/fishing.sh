#!/bin/bash

# Move to the data folder
cd ../data/

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 24 36 48)
# can use all treatments in the future, but for now we will just used Nutl because of memory issues on the computer.
# treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")

# Testing with Nutl only
treatments=("Nutl")
# Nutl does not have 36 hour time point so we should redefine this variable here
times=(0 8 12 16 24 48)


# salmon quant 
# -i salmon_index/Hsapiens_index_k17 
# -l A 
# -r trimmed/RPE_Nutl_t0_r1_tc.fastq 
# -p 8 
# --validateMappings 
# -o quantification/RPE_Nutl_t0_r1   

salmon index -t references/human_index_salmon/gentrome.fa -i references/human_index_salmon/Hsapiens_index_k17 --decoys references/human_index_salmon/decoys.txt -k 17 



# This command should run the salmon quant command for all treatments/ times.
for treatment in ${treatments[@]}; do
    for time in ${times[@]}; do
        time_folder="${treatment}_${time}"
        folder="organised/$treatment"

        echo "Trimming $time_folder"

        # folder for salmon outputs == aquarium obviously
        aquarium="organised/${treatment}/output_salmon"
        mkdir $aquarium

        input_trimmed_folder="organised/${treatment}/output_trimmed"
        input_trimmed_file="${input_trimmed_folder}/trimmed_fq_${time_folder}.fastq"

        salmon quant -i "references/human_index_salmon/Hsapiens_index_k17" \
        -l A \
        -r $input_trimmed_file \
        -p 8
        -o "$aquarium/salmon_quant_${time_folder}" \

        # Now we should delete the old file
        # rm $input_trimmed_file
    done
done
