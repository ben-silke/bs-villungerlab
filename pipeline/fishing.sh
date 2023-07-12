#!/bin/bash

echo "sleeping"
# sleep 10
echo "still sleeping"
# sleep 10000
echo "no longer sleeping"
# Move to the data folder
# cd ../data/
main_path=$1
cd $main_path

# Move to the data folder
# external_drive="/Volumes/bs_external/villunger"
# cd $external_drive

# Include all times, can exclude based upon experiement later.
# TODO: make these env variables so that you only need to update this once

times=(0 8 12 16 20 24 36 48)
replicates=("r1" "r2" "r3" "r4" "r5" "r6")
treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")

# SPECIFIC TREATMENT
# treatments=("ZM")
# ZM does not have 8 or 12 so can exclude
# times=(0 16 24 36 48)

treatments=("ZM" "Nutl" "Noc" "Etop" "DHCB")

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
# number=$(( $RANDOM % (30) ))

for treatment in ${treatments[@]}; do
    folder="organised/$treatment"
    aquarium="organised/${treatment}/output_salmon"
    mkdir $aquarium
    # sleep $number
    # number = $(($number + 1000))

    for time in ${times[@]}; do
        # number=$(($number + $time))
        # echo "Sleeping for $number"
        # sleep $number
        for replicate in ${replicates[@]}; do
            # Runs in parallel (phew)
            # (
            # echo "sleeping for $number"
            # sleep $number
            file="${treatment}_${time}_${replicate}"
            echo "Fishing in $file"

            if test -d "$aquarium/salmon_quant_${file}"; then
                echo "$aquarium/salmon_quant_${file} exists"
            else
                # folder for salmon outputs == aquarium obviously
                input_trimmed_folder="organised/${treatment}/output_trimmed"
                input_trimmed_file="${input_trimmed_folder}/trimmed_fq_${file}.fastq"

                salmon quant -i "references/human_index_salmon/Hsapiens_index_k17"  -o "$aquarium/salmon_quant_${file}" \
                -l A \
                -r $input_trimmed_file \
                -p 8

                sleep 100
            fi

            # ) &
            # Now we should delete the old file
            # rm $input_trimmed_file
        done
    done
done