'''
Missense variant miner:
find missense variants in ExAC
returns csv file.



'''
from __future__ import print_function

import sys, traceback, subprocess, gzip,  tarfile, os, signal
import pandas as pd
import csv

from pandas import DataFrame

import zipfile
from collections import OrderedDict
import collections

import string
import json
import pdb


VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

# Sets protein ID to search in dataframe
#ENSP = "ENSP00000003084" #CFTR
#ENSP = "ENSP00000269228" #NPC1
ENSP = "ENSP00000262304" #PKD1
ENSP = "ENSP00000262410" #MAPT
#ENSG = "ENSG00000001626" #CFTR
#ENSG = "ENSG00000141458" #NPC1
ENSG = "ENSG00000008710" #PKD1
ENSG = "ENSG00000186868" #ENSG
# ENSG = "ENSG00000186868" #MAPT
# ENSG = "ENSG00000272636" #Diagnostic - beginning of Chr17
#ENST = "ENST00000003084" CFTR
#ENST = 'ENST00000269228' #NPC1i
ENST = "ENST00000262304" #PKD1
ENST = "ENST00000262410" #MAPT


GENE = "MAPT"
FILENAME = "CFTR_PROV_extract.csv"
FILENAME1 = "CFTR_PROVEANScores.csv"
FILENAME2 = "MAPT_ExACScores.csv"
FILENAME3 = "MAPT_MutPredScores.csv"
FILENAME4 = "dbNSFP_output.csv"
FILENAME5 = "dbNSFP_extract.csv"
FILENAME6 = "MAPT_ExACScores.csv"
#UniProt = "P13569" CFTR
#UniProt = "O15118" #NPC1
UniProt = "P98161" #PKD1
UniProt = "P10636" #MAPT
#Chr = '18' NPC1
Chr = "16" #PKD1
Chr = "17"

# change directory to working with DAta
#os.chdir("../Data/")
cwd = os.getcwd()



def count_comments(filename):
    """Count comment lines (those that start with "#") in an optionally
	gzipped file.
	:param filename:  An optionally gzipped file.
	https://gist.github.com/slowkow/6215557
	"""
    comments = 0

    # fn_open = gzip.open if filename.endswith('.gz') else open
    # with fn_open(filename) as fh:
    for line in filename:
        if line.startswith('#'):
            comments += 1
            #print(line)
        else:

            break

    return comments


def parse(line):
    """Parse a single VCF line and return an OrderedDict.
	https://gist.github.com/slowkow/6215557
	"""
    result = OrderedDict()

    fields = line.rstrip().split('\t')

    # Read the values in the first seven columns.
    for i, col in enumerate(VCF_HEADER[:7]):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = fields[7].split(';')

    for i, info in enumerate(infos, 1):
        # info should be "key=value".
        try:
            key, value = info.split('=')
        # But sometimes it is just "value", so we'll make our own key.
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Set the value to None if there is no value.
        result[key] = _get_value(value)

    return result


def lines(Chr):
    """Open an optionally gzipped VCF file and generate an OrderedDict for
	each line.
	https://gist.github.com/slowkow/6215557
	"""
    # TODO: see if there is a way to first map chromosomes within file and keep this data in a temp file?
    #print('opening file')
    #fn_open = gzip.open #if filename.endswith('.gz') else open
    with open('/work/in/ExAC_data/ExAC.r0.3.1.sites.vep.vcf', 'rt') as fh, open(FILENAME6, 'w') as csvout:
        a = csv.writer(csvout)
        first_row = ('GENE_ID', 'TRANSCRIPT', 'TRANSCRIPT CHANGE', 'PROTEIN CHANGE', 'AA_POS', 'AA_CHANGE', 'MUTATION',
                     'ALLELE COUNT', 'ALLELE FREQUENCY')
        a.writerows([first_row])
        print ('opened file')
        print('looking for chromosome ', Chr)
        for filename in os.listdir(os.getcwd()):
		print(filename)
		cwd=os.getcwd()
		print cwd
	#items = list(range(1, 23))
        #l = len(items)
        FinalResults = OrderedDict()

        # printProgress(i, l, prefix = 'Progress:', suffix = 'Complete', barLength = 50)
        for good_line in find_good_lines(fh):
            # parse data line to find protein ID
            # dict item ['CQS'] contains ENSP IDs and mutation ID in tricode
            query = parse(good_line)['CSQ']

            # search for protein coding variants for POI
            if any(ENST in s for s in query):
                # get the rest of the data for matching lines
                p_good_line = parse(good_line)
                for k, e in enumerate(query):
                    for line in e.splitlines():
                        if ENSP in line and "missense_variant" in line:  # TODO: add deletions
                            AlleleCount = ()
                            AlleleFrequency = ()
                            #print(':::::')
                            # print(p_good_line['AF'])
                            #print(len(p_good_line['AC']))
                            if len(p_good_line['AC']) > 1:
                                #print('AC1', (p_good_line['AC']))
                                #print('AF1', (p_good_line['AF']))
                                try:
                                    AlleleCount = sum(list(map(int, p_good_line[
                                        'AC'])))  # list added to break object, if iterate dont need list(
                                    AlleleFrequency = sum(list(map(float, p_good_line['AF'])))
                                except ValueError:
                                    #print("error")
                                    AlleleCount = p_good_line['AC']
                                    AlleleFrequency = float(p_good_line['AF'])
                                    #print('AF2:', p_good_line['AF'])
                                    #print('AC2:', p_good_line['AC'])
                            else:
                                AlleleCount = int(p_good_line['AC'])
                                AlleleFrequency = float(p_good_line['AF'])

                            match = line.split( '|')  # GeneSym:match[3], geneID:match[4], ENST:match[6], ENST+change: match[10], ESNP_change:match [11], AApos:match[14], AA_change: match[15]
                            # print(match[3], match[6],match[10],match[11],match[14], match[15])
                            aa_change = match[15].split('/')
                            variant = [aa_change[0], match[14], aa_change[1]]
                            #print(variant)
                            # print(':::::::::::::::::::::::::::::::::::::::::::::::')
                            row_line = [match[3], match[6], match[10], match[11], match[14], match[15],
                                        ''.join(variant), AlleleCount, AlleleFrequency]
                            print(row_line)
                            a.writerows([row_line])
    # print(FinalResults)

    output_dict = json.loads(json.dumps(FinalResults))
    with open('ExacDUMP.txt', 'w') as outfile:
        json.dump(output_dict, outfile)


def find_good_lines(fh):
    #TODO: rewrite to account for Chr vs int variables and to reduce unnecessary searching.
    i='0'
    #print (fh)
    for line in fh:

        chrom = line[0:2]
        #print (chrom, Chr)
        if line.startswith('#'):
            continue
        if chrom == Chr:
            yield line
        else:
            if line[0:2] != i:
                #print (line[0:2])
                i = line[0:2]


            continue


        # if chrom > Chr:
        #     print ('end')
        #     break


def _get_value(value):
    """Interpret null values and return ``None``. Return a list if the value
	contains a comma.
	"""
    if not value or value in ['', '.', 'NA']:
        return None
    if ',' in value:
        return value.split(',')
    return value


def filterExACoutput(file):
    #print('reading file')
    df = pd.read_csv(file)
    col_names = []

    col_names = [i for i in string.printable[:len(df.columns)]]
    df.columns = [c.rstrip() for c in df.columns]
    #print(col_names)
    df.columns = [col_names]
    # df[1] = [col_names[1].astype(str)]
    # df[]


    #print('done')
    df = df[df['1'] == "missense_variant"]

    df = df[df['3'] == GENE]
    df = df[df['6'] == ENST]
    df.to_csv('Filtered_Exac_OUT.csv')

    #print(df.head())



def main():


    lines(Chr)
    #lines('ExAC.r0.3.1.sites.vep.vcf.gz', Chr)

    #filterExACoutput('Exac_parse_OUT.csv')




if __name__ == '__main__':
    main()
