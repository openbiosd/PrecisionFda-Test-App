'''
Missense variant miner:
find missense variants in ExAC, score them using PROVEAN, MutPred and dbNSFP
return csv file.



'''
from __future__ import print_function
import os
import sys
import pandas as pd
import csv

import zipfile




Chr = sys.argv[1]
ENSG = sys.argv[2]
OUT_FILE = "ExAC_OUT.csv"


FILENAME4 = "dbNSFP_output.csv"








def mine_dbNSFP(Chr, ENSG):
    # TODO: should be able to read file names and iterate through them
    chrfilesdict = {
        '1': 'dbNSFP3.2a_variant.chr1',
        '2': 'dbNSFP3.2a_variant.chr2',
        '3': 'dbNSFP3.2a_variant.chr3',
        '4': 'dbNSFP3.2a_variant.chr4',
        '5': 'dbNSFP3.2a_variant.chr5',
        '6': 'dbNSFP3.2a_variant.chr6',
        '7': 'dbNSFP3.2a_variant.chr7',
        '8': 'dbNSFP3.2a_variant.chr8',
        '9': 'dbNSFP3.2a_variant.chr9',
        '10': 'dbNSFP3.2a_variant.chr10',
        '11': 'dbNSFP3.2a_variant.chr11',
        '12': 'dbNSFP3.2a_variant.chr12',
        '13': 'dbNSFP3.2a_variant.chr13',
        '14': 'dbNSFP3.2a_variant.chr14',
        '15': 'dbNSFP3.2a_variant.chr15',
        '16': 'dbNSFP3.2a_variant.chr16',
        '17': 'dbNSFP3.2a_variant.chr17',
        '18': 'dbNSFP3.2a_variant.chr18',
        '19': 'dbNSFP3.2a_variant.chr19',
        '20': 'dbNSFP3.2a_variant.chr20',
        '21': 'dbNSFP3.2a_variant.chr21',
        'X': 'dbNSFP3.2a_variant.chrX',
        'M': 'dbNSFP3.2a_variant.chrM',
    }
    # read from tsv.gz file/work/in/ExAC_data/ExAC.r0.3.1.sites.vep.vcf'
    with zipfile.ZipFile('/work/in/dbnsfp/dbNSFPv3.2a.zip', 'r') as tsvin, open(OUT_FILE, 'wt') as csvout:
        tp = pd.read_csv(tsvin.open(chrfilesdict[Chr]), delimiter='\t', quoting=csv.QUOTE_NONE, iterator=True,
                         dtype=object, chunksize = 100000)  # header = None)

        writer = csv.writer(csvout)
        row1 = ['#chr	pos(1-based)	ref	alt	aaref	aaalt	rs_dbSNP146	hg19_chr	hg19_pos(1-based)	hg18_chr	hg18_pos(1-based)	genename	cds_strand	refcodon	codonpos	codon_degeneracy	Ancestral_allele	AltaiNeandertal	Denisova	Ensembl_geneid	Ensembl_transcriptid	Ensembl_proteinid	aapos	SIFT_score	SIFT_converted_rankscore	SIFT_pred	Uniprot_acc_Polyphen2	Uniprot_id_Polyphen2	Uniprot_aapos_Polyphen2	Polyphen2_HDIV_score	Polyphen2_HDIV_rankscore	Polyphen2_HDIV_pred	Polyphen2_HVAR_score	Polyphen2_HVAR_rankscore	Polyphen2_HVAR_pred	LRT_score	LRT_converted_rankscore	LRT_pred	LRT_Omega	MutationTaster_score	MutationTaster_converted_rankscore	MutationTaster_pred	MutationTaster_model	MutationTaster_AAE	MutationAssessor_UniprotID	MutationAssessor_variant	MutationAssessor_score	MutationAssessor_score_rankscore	MutationAssessor_pred	FATHMM_score	FATHMM_converted_rankscore	FATHMM_pred	PROVEAN_score	PROVEAN_converted_rankscore	PROVEAN_pred	Transcript_id_VEST3	Transcript_var_VEST3	VEST3_score	VEST3_rankscore	MetaSVM_score	MetaSVM_rankscore	MetaSVM_pred	MetaLR_score	MetaLR_rankscore	MetaLR_pred	Reliability_index	CADD_raw	CADD_raw_rankscore	CADD_phred	DANN_score	DANN_rankscore	fathmm-MKL_coding_score	fathmm-MKL_coding_rankscore	fathmm-MKL_coding_pred	fathmm-MKL_coding_group	Eigen-raw	Eigen-phred	Eigen-raw_rankscore	Eigen-PC-raw	Eigen-PC-raw_rankscore	GenoCanyon_score	GenoCanyon_score_rankscore	integrated_fitCons_score	integrated_fitCons_score_rankscore	integrated_confidence_value	GM12878_fitCons_score	GM12878_fitCons_score_rankscore	GM12878_confidence_value	H1-hESC_fitCons_score	H1-hESC_fitCons_score_rankscore	H1-hESC_confidence_value	HUVEC_fitCons_score	HUVEC_fitCons_score_rankscore	HUVEC_confidence_value	GERP++_NR	GERP++_RS	GERP++_RS_rankscore	phyloP100way_vertebrate	phyloP100way_vertebrate_rankscore	phyloP20way_mammalian	phyloP20way_mammalian_rankscore	phastCons100way_vertebrate	phastCons100way_vertebrate_rankscore	phastCons20way_mammalian	phastCons20way_mammalian_rankscore	SiPhy_29way_pi	SiPhy_29way_logOdds	SiPhy_29way_logOdds_rankscore	1000Gp3_AC	1000Gp3_AF	1000Gp3_AFR_AC	1000Gp3_AFR_AF	1000Gp3_EUR_AC	1000Gp3_EUR_AF	1000Gp3_AMR_AC	1000Gp3_AMR_AF	1000Gp3_EAS_AC	1000Gp3_EAS_AF	1000Gp3_SAS_AC	1000Gp3_SAS_AF	TWINSUK_AC	TWINSUK_AF	ALSPAC_AC	ALSPAC_AF	ESP6500_AA_AC	ESP6500_AA_AF	ESP6500_EA_AC	ESP6500_EA_AF	ExAC_AC	ExAC_AF	ExAC_Adj_AC	ExAC_Adj_AF	ExAC_AFR_AC	ExAC_AFR_AF	ExAC_AMR_AC	ExAC_AMR_AF	ExAC_EAS_AC	ExAC_EAS_AF	ExAC_FIN_AC	ExAC_FIN_AF	ExAC_NFE_AC	ExAC_NFE_AF	ExAC_SAS_AC	ExAC_SAS_AF	ExAC_nonTCGA_AC	ExAC_nonTCGA_AF	ExAC_nonTCGA_Adj_AC	ExAC_nonTCGA_Adj_AF	ExAC_nonTCGA_AFR_AC	ExAC_nonTCGA_AFR_AF	ExAC_nonTCGA_AMR_AC	ExAC_nonTCGA_AMR_AF	ExAC_nonTCGA_EAS_AC	ExAC_nonTCGA_EAS_AF	ExAC_nonTCGA_FIN_AC	ExAC_nonTCGA_FIN_AF	ExAC_nonTCGA_NFE_AC	ExAC_nonTCGA_NFE_AF	ExAC_nonTCGA_SAS_AC	ExAC_nonTCGA_SAS_AF	ExAC_nonpsych_AC	ExAC_nonpsych_AF	ExAC_nonpsych_Adj_AC	ExAC_nonpsych_Adj_AF	ExAC_nonpsych_AFR_AC	ExAC_nonpsych_AFR_AF	ExAC_nonpsych_AMR_AC	ExAC_nonpsych_AMR_AF	ExAC_nonpsych_EAS_AC	ExAC_nonpsych_EAS_AF	ExAC_nonpsych_FIN_AC	ExAC_nonpsych_FIN_AF	ExAC_nonpsych_NFE_AC	ExAC_nonpsych_NFE_AF	ExAC_nonpsych_SAS_AC	ExAC_nonpsych_SAS_AF	clinvar_rs	clinvar_clnsig	clinvar_trait	clinvar_golden_stars	Interpro_domain	GTEx_V6_gene, GTEx_V6_tissue']
        # for i in range(1):
        #     row1 = next(tp)
        #     print(row1)
        #     print('found ', len(row1), 'rows')
        #     writer.writerows([row1[1:len(row1)]])
        #     i = +1
        # print(chrfilesdict[Chr])
        # tp = pd.read_csv(tsvin.open(chrfilesdict[Chr]), delimiter='\t', quoting=csv.QUOTE_NONE ,iterator=True,
        #                  dtype=object, chunksize=40)  # header = None)
        print ("Searching Chromosome:", Chr, "ENSG:", ENSG)
        variants = 0

        #TODO: see if there is a way to get the first row of every chunk, to speed up search
        for chunk in tp:
            row = next(chunk.itertuples())
            count = row[20]

            #print (row[1:len(row)])
            if count == ENSG:
                variants += 1
                print(variants, ' writing', ENSG, 'variant scores ', row[0])
                writer.writerows([row[1:len(row)]])
        print(variants)


def db_NSFP_iterate(fh):
    #TODO: rewrite to account for Chr vs int variables and to reduce unnecessary searching.
    i='0'
    found = 0

    for line in fh:
        print (line)

        ENSG_match = line[20]

        if ENSG_match == ENSG:
            print (line)
            yield line
            found = 1
        else:
            if ENSG_match != ENSG and found ==1:
                print ("no more to be found")
                break
            continue

def cleanup_dbNSFP_extract(file):

    with open(file, 'rt') as tsvin, open(OUT_FILE, 'wt') as csvout:
        dict = []
        df = pd.read_csv(tsvin, delimiter=',', encoding="utf-8-sig")

        for row in enumerate(df['FATHMM_score']):

            FAS = row[1]

            print (FAS)
            if FAS[0] == FAS[1]:
                id = row[0]
                print(df.iloc[[id]])


def extract_dbNSFP(file):

    """ There are rows with multiple comma separated values in them
        this function works to convert such rows into rows which have
        one value per line
        Problems is that it will not run unless rows have exact number of comma sep. values
        so next step is to subset these rows and then merge with the rest of the dataframe
    """
    with open(file, 'rt') as tsvin, open(OUT_FILE, 'wt') as csvout:
        dict = []
        df = pd.read_csv(tsvin, delimiter=',', encoding="utf-8-sig")

        #print (df.head(1))
        #print (df['Ensembl_transcriptid'])
        #column_names = ['Ensembl_transcriptid', 'Ensembl_proteinid',  'MutationTaster_score','MutationTaster_pred', 'MutationTaster_AAE', 'FATHMM_score', 'FATHMM_pred']
        #print(column_names)
        df = df[['Ensembl_transcriptid', 'Ensembl_proteinid',  'MutationTaster_score','MutationTaster_pred', 'MutationTaster_AAE', 'FATHMM_score', 'FATHMM_pred']]

        for col in ['Ensembl_transcriptid', 'Ensembl_proteinid',  'MutationTaster_score','MutationTaster_pred', 'MutationTaster_AAE', 'FATHMM_score', 'FATHMM_pred']:
            df [col] = df[col].str.split(';')



            #print (len(row['Ensembl_proteinid']))'Ensembl_proteinid'

'Ensembl_proteinid'
        #
        # #print (df['Ensembl_transcriptid'])
        # #df.to_csv(csvout)
        # i =  df['Ensembl_transcriptid'].map(len)
        # j = np.repeat(np.arange(len(df)),i)
        # k = np.concatenate(list(map(np.arange, i)))
        # df = df.iloc[j]
        # print (df['Ensembl_transcriptid'])
        # for col in ['Ensembl_transcriptid', 'Ensembl_proteinid',  'MutationTaster_score','MutationTaster_pred', 'MutationTaster_AAE', 'FATHMM_score', 'FATHMM_pred']:
        #     df [col] = list(map(lambda xs, i: xs[i], df[col], k))
        # df.to_csv(csvout)




def main():

 mine_dbNSFP(Chr, ENSG)
# extract_dbNSFP(FILENAME4)



if __name__ == '__main__':
    main()