#!/bin/bash

# folder names can be created better but this is just to test

main_path=$1
cd $main_path

# Move to the data folder
# cd ../data/raw/

# Move to the data folder
# external_drive="/Volumes/bs_external/villunger"
# cd "$external_drive/raw"
# mkdir ../organised

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 20 24 36 48)
replicates=("r1" "r2" "r3")

# can use all treatments in the future, but for now we will just used Nutl because of space issues on the computer.
# treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")
treatments=("ZM" "Nutl" "Noc" "Etop" "DHCB")

# treatments=("ZM")
# ZM does not have 8 or 12 so can exclude
# times=(0 16 24 36 48)

for treatment in ${treatments[@]}; do
    echo "treatment = $treatment"
    mkdir ../organised/$treatment
    for time in ${times[@]}; do
        echo "time=$time"
        time_folder="${treatment}_${time}"
        echo "processing $time_folder"
        mkdir ../organised/$treatment/$time_folder
        
        for replicate in ${replicates[@]}; do
            replicate_folder="../organised/$treatment/$time_folder/$replicate"
            mkdir $replicate_folder
            mkdir "$replicate_folder/experiment_one"
            mkdir "$replicate_folder/experiment_two"
            # Code block to find all treatment and times
            # This should be reasonably non-specific and not cause issues if the directory structure is slightly different;
            # provided that .*samples is within the folder name.
            find . -type d -name "*samples" | while read -r directory; do
                # Loop through each file matching the pattern
                echo "processing $directory"
                file_expression="#${treatment}_t${time}_${replicate}"
                for file in $directory/*$file_expression*.bam; do
                    # Move the file into the new folder
                    # echo $file
                    substring="${file#*/}"  # Removes everything before the first slash
                    experiment_folder="${substring%%/*}"
                    # echo "experiment_folder $experiment_folder"
                    file_name="${file##*/}"
                    # echo "file_name = $file_name"
                    # echo "$replicate_folder/$experiment_folder/$file_name"
                    cp "$file" "$replicate_folder/$experiment_folder/$file_name"
                    # Could also use mv, depending on how you feel about deleting files; cp for now and the files in that subsequent directory will be removed after merge
                done
            done
        done
    done
done