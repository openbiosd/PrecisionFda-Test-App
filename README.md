# PrecisionFda-Test-App

Parser for ExAC data which describes genetic variants in >60,000 individuals. 
For more information on ExAC read here:http://exac.broadinstitute.org/about
Please read the McArthur blog for even more information: https://macarthurlab.org/blog/
Also listen to this presentation by Daniel McArthur about human variation and ExAC https://youtu.be/bM977g4hOz0

usage: exac_parser.py Chr, ENST, OUT_FILE.vcf

example: exac_parser.py 17, ENST00000262410, MAPT_OUT.vcf


exac_parser.py returns a csv file with missense variants for specified gene transcripts:
gene ID, transcript ID, variant position (amino acid), reference (REF) amino acid, variant amino acid (ALT), allele frequncy (AF) and Allele Count (AC) 
but this can be extended to much more, such as population metrics, splicing variants, SNPs etc etc etc. 

function <i>search_exac(Chr, ENST, OUT_FILE)</i>: 
requires a local copy of current ExAC release (r0.3.1) on your machine in ~/Data/ to run. 
Download from here (4GB): ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz 
for readme on Exac file look in ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/



