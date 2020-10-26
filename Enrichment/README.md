binomial_run.py uses the R script binomial.R and one file and one directory as input: 

   INPUT_1_portal_orthogroup_counts.tsv = with columns (portal group total core  group-specific acc) 
   INPUT_go_freq_examples = directory of *.tsv files, split from the frequencies tools output by portal

and will output: 

   OUTPUT.TSV = with adjuested pvalues for each term in the portal in columns (portal  term  obs1 total1 obs2 total2 p_val adj_p_val)
