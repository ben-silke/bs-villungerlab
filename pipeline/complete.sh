
# Steps:

# 1. move files - moves files from the raw directory to grouped directories
# 2. merge_bam_files.sh - merges the two lanes into a single bam file
# 3. convert_bam_to_fastq.sh - converts the bam files to fastq format.
# 4. trim_reads.sh - trims the fastq file formats
# 5. fishing/ STAR Gazing, - this step is dependent upon the analysis which you would like to run; fishing uses salmon, while STAR gazing uses STAR


./move_files.sh /Volumes/bs_external/villunger
./merge_bam_files.sh /Volumes/bs_external/villunger
./convert_bam_to_fastq.sh /Volumes/bs_external/villunger
./trim_reads.sh /Volumes/bs_external/villunger
./fishing.sh /Volumes/bs_external/villunger