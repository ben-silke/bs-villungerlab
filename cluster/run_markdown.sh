#!/bin/bash

Rscript ../r/src/salmon_DESeq2/local_files/DHCB_r1to6_create_data.R
Rscript ../r/src/salmon_DESeq2/local_files/ZM_r1to6_create_data.R
Rscript ../r/src/salmon_DESeq2/local_files/Etop_r1to6_create_data.R
Rscript ../r/src/salmon_DESeq2/local_files/Nutl_r1to6_create_data.R
Rscript ../r/src/salmon_DESeq2/local_files/Noc_r1to6_create_data.R

sleep 1000

Rscript ../r/src/salmon_DESeq2/local_files/knit_all_files.R