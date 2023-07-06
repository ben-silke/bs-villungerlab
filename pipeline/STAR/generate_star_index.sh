
# Index generation
# STAR \
# --runThreadN 2 \
# --runMode genomeGenerate \
# --genomeDir ../../data/references/STAR \
# --genomeFastaFiles ../../data/references/human_index_salmon/gentrome.fa

# Running mapping
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir data/references/STAR --genomeFastaFiles data/references/human_index_salmon/gentrome.fa