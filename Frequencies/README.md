get_orthogroup_annotation_freq.py uses two files as input:
   INPUT_2_portal_orthogroup_group_id.tsv = with columns (portal  group   og      id)
   INPUT_functional_annotation_GO.tsv = with concatentated annotation files for all species

and will output:
   OUTPUT.TSV = with occurances by orthogroup in columns (portal term  total core group-specific(mixed)  acc)
   
Enrichment analysis required that the OUTPUT.TSV file be broken up by portal using large_file_split.py
 
