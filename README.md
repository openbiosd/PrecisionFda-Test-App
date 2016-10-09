# PrecisionFda-Test-App

Parser for ExAC data which describes genetic variants in >60,000 individuals. 
For more information on ExAC read here:http://exac.broadinstitute.org/about
Please read the McArthur blog for even more information: https://macarthurlab.org/blog/

exac_parser.py is currently configured to mine missense variants for specific gene transcripts and proteins for:
gene ID, transcript ID, variant position (amino acid), reference (REF) amino acid, variant amino acid (ALT), allele frequncy (AF) and Allele Count (AC) 
but this can be extended to much more, such as population metrics, splicing variants, SNPs etc etc etc. 

#<i>lines(ENST,ENSP,Chr)</i>: 
takes in ENST, ENSP and Chr (transcript id, protein id and chromosome #)  and returns a specified .csv file with missense variants in the queried protein. Requires a local copy of current ExAC release (r0.3.1) on your machine in ~/Data/ to run. 
Download from here (4GB): ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz 
for readme on Exac file look in ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/

other functions which are working (badly) and are commented out in exac_parser.py are: 

#<i>mine_PROVEAN:</i>
mines pre-computed consequence scores for amino acid substitutions for all known proteins
http://provean.jcvi.org/index.php

#<i>mine_dbNSFP:</i>
mines an aggregate of pre-compute consequence scores for all known protein
https://sites.google.com/site/jpopgen/dbNSFP

#<i>mine_MutPred:</i> mines pre-computed consequence scores and biochemical hypothesis for all known proteins
http://mutpred.mutdb.org/about.html

These are optional feautures for this project, but it would be great to join the output of all of these into one dataframe in python.
