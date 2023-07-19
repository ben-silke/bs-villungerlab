import os

with open('edited_gtf_annotation_file.gtf .gtf', 'r') as fin, open('output_updated_annotations.gtf', 'w') as fout:
    for line in fin:
        fields = line.strip().split(';')
        transcript_id = [f.strip().split(' ')[1].replace('"', '') for f in fields if 'transcript_id' in f][0]
        transcript_version = [f.strip().split(' ')[1].replace('"', '') for f in fields if 'transcript_version' in f][0]
        new_line = line.strip() + '; transcript_vid "{}.{}";\n'.format(transcript_id, transcript_version)
        print(f'{new_line=}')
        fout.write(new_line)



# import os

# with open('edited_gtf_annotation_file.gtf', 'r') as fin, open('output_updated_annotations.gtf', 'w') as fout:
#     for line in fin:
#         print(line)
#         fields = line.strip().split(';')
#         print(fields)
#         if len(fields) > 4 and 'transcript_id' in line:
#             transcript_id = [f.strip().split(' ')[1].replace('"', '') for f in fields if 'transcript_id' in f][0]
#             transcript_version = [f.strip().split(' ')[1].replace('"', '') for f in fields if 'transcript_version' in f][0]
#             new_line = line.strip() + ' transcript_vid "{}.{}";\n'.format(transcript_id, transcript_version)
#             print(new_line)
#             fout.write(new_line)