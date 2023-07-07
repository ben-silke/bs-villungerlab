#!/bin/bash

# To retain ability to run on different files.
BBMAP_PATH=$1


# bbduk.sh 
# in=raw_data/RPE1_Nutl_t0_r1.fastq 
# out=trimmed/RPE_Nutl_t0_r1_tc.fastq 
# ref=/home/dario/bbmap/resources/polyA.fa.gz,/home/dario/bbmap/resources/truseq.fa.gz 
# k=13 
# ktrim=r 
# useshortkmers=t 
# mink=5 
# qtrim=t 
# trimq=10 
# minlength=20


# Move to the data folder
# cd ../data/
external_drive="/Volumes/bs_external/villunger"
cd $external_drive

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 24 36 48)
treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")
replicates=("r1" "r2" "r3" "r4" "r5" "r6")

## SPECIFIC TREATMENT
# treatments=("Nutl")
# # Nutl does not have 36 hour time point so we should redefine this variable here
# times=(0 8 12 16 24 48)

# SPECIFIC TREATMENT
treatments=("ZM")
# ZM does not have 8 or 12 so can exclude
times=(0 16 24 36 48)

# This command should run the bbduk decomtamination for all treatments/ times.
for treatment in ${treatments[@]}; do
    for time in ${times[@]}; do
        for replicate in ${replicates[@]}; do

            file="${treatment}_${time}_${replicate}"

            echo "Trimming $file"
            output_fastq_folder="organised/${treatment}/output_merged_fq"
            output_trimmed_folder="organised/${treatment}/output_trimmed"
            mkdir $output_trimmed_folder
            
            # Command to run bbduk.sh with necessary arguments.
            # See https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/ for arguments
            $BBMAP_PATH/bbduk.sh \
            in="$output_fastq_folder/output_merged_fq_${file}.fastq" \
            out="${output_trimmed_folder}/trimmed_fq_${file}.fastq" \
            stats="../project-plan/validation/${treatment}/stats_${file}.txt" \
            # This is just removing a polyA tail and the primer, rRNA is not being removed here
            ref="$BBMAP_PATH/resources/polyA.fa.gz,$BBMAP_PATH/resources/truseq.fa.gz" \
            k=13 \
            ktrim=r \
            useshortkmers=t \
            mink=5 \
            qtrim=t \
            trimq=10 \
            minlength=20 \

            # Now we should delete the old file
            rm "$output_fastq_folder/output_merged_fq_${file}.fastq"
        done
    done
done
