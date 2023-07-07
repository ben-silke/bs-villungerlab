#!/bin/bash

# Move to the data folder
# cd ../data/

main_path=$1
cd $main_path
ls
pwd
# exit
# Move to the data folder
# external_drive="/Volumes/bs_external/villunger"
# cd $external_drive


# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 20 24 36 48)
replicates=("r1" "r2" "r3")
experiment_folders=("experiment_one" "experiment_two")

# can use all treatments in the future, but for now we will just used Nutl because of memory issues on the computer.
treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")
treatments=("ZM" "Nutl" "Noc" "Etop" "DHCB")


# SPECIFIC TREATMENT
# treatments=("ZM")
# ZM does not have 8 or 12 so can exclude
# times=(0 16 24 36 48)

# This command should run the samtools merge command for all treatments/ times.
for treatment in ${treatments[@]}; do
    output_merged_folder="organised/${treatment}/output_merged"
    mkdir $output_merged_folder
    for time in ${times[@]}; do
        time_folder="${treatment}_${time}"

        for replicate in ${replicates[@]}; do
            for experiment_folder in ${experiment_folders[@]}; do

                folder="organised/$treatment/$time_folder/$replicate/$experiment_folder"
                echo "Merging $folder"
                if [ "$experiment_folder" == "experiment_one" ]; then
                    file_name="$replicate"
                else
                    if [ "$replicate" == "r1" ]; then
                        file_name="r4"
                    elif [ "$replicate" == "r2" ]; then
                        file_name="r5"
                    else 
                        file_name="r6"
                    fi
                fi
                
                echo "Merging contents of $folder into $output_merged_folder/output_merged_${time_folder}_${file_name}"
                file_name="$output_merged_folder/output_merged_${time_folder}_${file_name}.bam"
                # NOTE: We may need to look at sorting these files before using this command.
                samtools merge -o $file_name $folder/*.bam

                # To remove the files
                # These files are copied from data/raw so can easily be recopied.
                # rm $folder/*.bam
                # rm $folder
            done
        done
    done
done