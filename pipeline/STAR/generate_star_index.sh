
# Index generation
# STAR \
# --runThreadN 2 \
# --runMode genomeGenerate \
# --genomeDir ../../data/references/STAR \
# --genomeFastaFiles ../../data/references/human_index_salmon/gentrome.fa

# /Volumes/bs_external/villunger/references/STAR/ensembl_genome/combined_reference.fa

# Running mapping
external_drive="/Volumes/bs_external/villunger/references"
echo "$external_drive"

# Download the ensembl genome
# rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-109/fasta/homo_sapiens/dna "$external_drive/STAR/ensembl_genome"
# Download the annotation
# rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-109/gff3/homo_sapiens "$external_drive/STAR/ensemble_annotation"

# merge the cDNA and the ncRNA


# This requires 32 GB ram and therefore needs to be run on the cluster
STAR \
--runThreadN 2 \
--runMode genomeGenerate \
--genomeDir "$external_drive/STAR/star_index" \
--genomeFastaFiles "$external_drive/STAR/ensembl_genome/combined_reference.fa" \
--sjdbGTFfile "$external_drive/STAR/ensemble_annotation" \
