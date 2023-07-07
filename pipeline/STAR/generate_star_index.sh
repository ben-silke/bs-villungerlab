
# Index generation
# STAR \
# --runThreadN 2 \
# --runMode genomeGenerate \
# --genomeDir ../../data/references/STAR \
# --genomeFastaFiles ../../data/references/human_index_salmon/gentrome.fa

# Running mapping
external_drive="/Volumes/bs_external/villunger/references"
echo "$external_drive"

# Download the ensembl genome
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-109/fasta/homo_sapiens/dna "$external_drive/STAR/ensembl_genome"
# Download the annotation
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-109/gff3/homo_sapiens "$external_drive/STAR/ensemble_annotation"


STAR \
--runThreadN 2 \
--runMode genomeGenerate \
--genomeDir "$external_drive/STAR/star_index" \
--genomeFastaFiles "$external_drive/STAR/ensembl_genome/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
--sjdbGTFfile "$external_drive/STAR/ensemble_annotation" \ # path to file which contains reference
date
