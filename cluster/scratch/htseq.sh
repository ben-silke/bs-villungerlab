 1042  htseq -i "transcript_vid" ZM_0_r1Aligned.out.sam ../references/output_updated_annotations.gtf
 1043  htseq-count -i "transcript_vid" ZM_0_r1Aligned.out.sam ../references/output_updated_annotations.gtf
 1044  htseq-count -i "transcript_vid" ZM_0_r1Aligned.out.sam ../../references/output_updated_annotations.gtf
 1045  ls
 1046  htseq-count --nonunique all -i "transcript_vid" ZM_0_r1Aligned.out.sam ../../references/output_updated_annotations.gtf