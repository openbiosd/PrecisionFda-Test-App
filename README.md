# PrecisionFda-Test-App

Parser for ExAC data which describes genetic variants in >60,000 individuals. 
For more information on ExAC read here:http://exac.broadinstitute.org/about
Please read the McArthur blog for even more information: https://macarthurlab.org/blog/

exac_parser.py is currently configured to mine missense variants for specific gene transcripts and proteins for:
gene ID, transcript ID, variant position (amino acid), reference (REF) amino acid, variant amino acid (ALT), allele frequncy (AF) and Allele Count (AC) 
but this can be extended to much more, such as population metrics, splicing variants, SNPs etc etc etc. 

parser_exac.py exac mining function is currently named lines. 

other functions which are currently working (badly) in exac_parser.py are: 

mine_PROVEAN: mines pre-computed consequence scores for amino acid substitutions for all known proteins
mine_dbNSFP: mines an aggregate of pre-compute consequence scores for all known protein
mine_MutPred: mines pre-computed consequence scores and biochemical hypothesis for all known proteins

These are optional feautures for this project, but it would be great to join the output of all of these into one dataframe in python.
