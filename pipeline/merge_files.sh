#!/bin/bash

# folder names can be created better but this is just to test


# Move to the data folder
cd ../data/raw/
mkdir ../organised

# Include all times, can exclude based upon experiement later.
times=(0 8 12 16 24 36 48)
# can use all treatments in the future, but for now we will just used Nutl because of space issues on the computer.
# treatments=("ZM" "Nutl" "Noc" "Nalm6_ZM" "Etop" "DHCB")
# Question -> Nalm6_ZM?
treatments=("Nutl")

for treatment in ${treatments[@]}; do
    echo "treatment = $treatment"
    mkdir ../organised/$treatment
    for time in ${times[@]}; do
        echo "time=$time"
        time_folder="${treatment}_${time}"
        echo $time_folder
        mkdir ../organised/$treatment/$time_folder
        
        
        # Code block to find all treatment and times
        # This should be reasonably non-specific and not cause issues if the directory structure is slightly different;
        # provided that .*samples is within the folder name.
        find . -type d -name "*samples" | while read -r directory; do
            # Loop through each file matching the pattern
            echo "processing $directory"
            file_expression="${treatment}_t${time}"
            for file in $directory/*$file_expression*.bam; do
                # Move the file into the new folder
                file_name="${file##*/}"
                echo "file_name = $file_name"
                cp "$file" "../organised/$treatment/$time_folder/$file_name"
                # Could also use mv, depending on how you feel about deleting files; cp for now and the files in that subsequent directory will be removed after merge
            done
        done
    done
done