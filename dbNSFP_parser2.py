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








def mine_dbNSFP(Chr, ENSG, OUT_FILE):
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
    with zipfile.ZipFile('/work/in/dbnsfp/dbNSFPv3.2a.zip', 'r') as tsvin, open(OUT_FILE, 'w') as csvout:
        tp = pd.read_csv(tsvin.open(chrfilesdict[Chr]), delimiter='\t', quoting=csv.QUOTE_NONE, iterator=True,
                         dtype=object, chunksize = 1000)  # header = None)

        writer = csv.writer(csvout)
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

 mine_dbNSFP(Chr, ENSG, OUT_FILE)
# extract_dbNSFP(FILENAME4)



if __name__ == '__main__':
    main()