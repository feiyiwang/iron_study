import pandas as pd
import numpy as np
import tqdm

# All the necessary paths
# path to data in red_library
freq_path = '/finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9.afreq'
event_path = '/finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_detailed_longitudinal_1.0.txt.gz'
pca_path = '/finngen/library-red/finngen_R9/pca_1.0/data/finngen_R9.eigenvec.txt'
sex_path = '/finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_minimum_1.0.txt.gz'
# path to uploaded data
snp_path = '/finngen/green/FeiyiWang/all_variants_and_proxies.csv'
phecode_path = '/finngen/green/FeiyiWang/phecode_map.csv'


# Obtain more accurate allele frequencies
# load data
freq = pd.read_csv(freq_path, sep='\t')
snp = pd.read_csv(snp_path)
# left join snp and freq
snp = snp.merge(freq[['ID','ALT_FREQS']], left_on='sandbox_format', right_on='ID')
snp = snp.rename(columns={'ALT_FREQS':'sandbox_af'})
# drop unnecessary cols
snp = snp.drop(columns=['Unnamed: 0', 'ID'])
# save the updated snp dataframe
snp.to_csv('all_variants_and_proxies.csv', index=None)


# Build a matrix of genetic dosages
# 1. save a list of SNP positions for obtaining genetic data at individual level via plink
snp_sub = snp[(~snp.sandbox_format.isna())&(snp.match == True)&(~snp.sandbox_format.duplicated())]
snp_list = snp_sub.sandbox_format.tolist()
snp_strings = ''
for i in snp_list:
    snp_strings += i+','

# 2. run plink to obatin a .ped file directly in notebook
#!plink --bfile /finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9 --snps $snp_strings --recode --out selected_snp
ped_head = {
    0:'family_id',
    1:'finngen_id',
    2:'father_id',
    3:'mother_id',
    4:'sex',
    5:'phenotype'
}
# In .ped file, the first 6 cols are listed as ped_head
# and then each two cols in the rest are a pair of alleles tied to a specific SNP
ped = pd.read_csv('selected_snp.ped', sep=' ', header=None).rename(columns=ped_head)

# 3. build an exposure matrix
# get a list of reference alleles
snp_ref = snp_sub.finn_ref.tolist()
# set up a dataframe with ids for the dosage matrix
exposure_matrix = ped[['finngen_id']]
# add sex to the dataframe
demo = pd.read_csv(sex_path, sep='\t')
exposure_matrix = exposure_matrix.merge(demo[['FINNGENID', 'SEX']], 'left', left_on='finngen_id', right_on='FINNGENID')
exposure_matrix = exposure_matrix.rename(columns={'SEX':'sex'})
exposure_matrix = exposure_matrix[['finngen_id', 'sex']]
# add dosages to the dataframe
for i in tqdm.tqdm(range(len(snp_list))):
    allele1 = np.select([(ped[6+2*i] == snp_ref[i]), (ped[6+2*i] != snp_ref[i])], [1, 0])
    allele2 = np.select([(ped[6+2*i+1] == snp_ref[i]), (ped[6+2*i+1] != snp_ref[i])], [1, 0])
    dosage = allele1 + allele2
    exposure_matrix[snp_list[i]] = dosage
# get read for the confounding matrix
confounding = exposure_matrix[['finngen_id', 'sex']]
# shape of exposure_matrix is 392649 by 2459
# 392649 rows -> 392649 individuals collected in R9
# 2459 cols -> 1 finngen_id + 2458 dosages
exposure_matrix = exposure_matrix.drop(columns=['sex'])
exposure_matrix.to_csv('exposure_matrix.csv', index=None)
# delete dataframes to release some memory
del ped, freq, snp, demo


# Build a matrix of 10 principal components of ancestry
pca = pd.read(pca_path)
confounding = confounding.merge(pca.iloc[:, 1:12], 'left', left_on='finngen_id', right_on='IID')
# shape of exposure_matrix is 392649 by 12
# 392649 rows -> 392649 individuals collected in R9
# 12 cols -> 1 finngen_id + 1 sex + 10 principal components
confounding = confounding.drop(columns=['IID'])
confounding.to_csv('confounding.csv', index=None)
# delete dataframes to release some memory
del pca


# Convert ICD codes to a list of phecodes
# 1. load data
events = pd.read_csv(event_path, sep='\t')
phecode_map = pd.read_csv(phecode_path)

# 2. clean the event data
# keep only icd9 and icd10
events = events[events.ICDVER.isin(['9','10'])][['FINNGENID', 'SOURCE', 'EVENT_AGE', 'CODE1']]
# remove n/a value
events = events[~events.CODE1.isna()]
# remove duplicates
events = events[~events.duplicated(subset=['FINNGENID', 'CODE1'])]
# split datasets by icd version
events_icd9 = events[events.ICDVER == '9']
events_icd10 = events[events.ICDVER == '10']
icd9 = phecode_map[phecode_map.vocabulary_id == 'ICD9CM'][['code', 'phecode']]
icd10 = phecode_map[phecode_map.vocabulary_id == 'ICD10CM'][['code', 'phecode']]
# clean event data coded by icd9
# 1) fix code with '}'
events_icd9.loc[34308875, 'CODE1'] = '4059'
# 2) remove codes with '_'
events_icd9 = events_icd9[events_icd9.CODE1 != '____']
# 3) fix codes with ':'
events_icd9.loc[74832158, 'CODE1'] = '3018'
events_icd9.loc[93924451, 'CODE1'] = '2152'
# 4) fix codes with '-'
to_fix = events_icd9[events_icd9.CODE1.str[-1] == '-'].CODE1.str[:-1]
for i, item in to_fix.iteritems():
    events_icd9.loc[i, 'CODE1'] = item
# 5) remove codes with '*'
events_icd9 = events_icd9[events_icd9.CODE1.str[-1] != '*']
# 6) remove code with '"'
events_icd9 = events_icd9.drop(76654182)
# 7) replace codes ending with lower case with upper case
to_fix = events_icd9[events_icd9.CODE1.str[-1].isin(['a', 'b', 'x'])].CODE1
for i, item in to_fix.iteritems():
    events_icd9.loc[i, 'CODE1'] = item.upper()

# 3. preprocess the datasets
def getSign(data):
    sign = []
    for i in tqdm.tqdm(data.icd):
        if len(i) == 5:
            df = data[data.icd.str.startswith(i[:-1])]
            if len(df) == 1:
                sign.append(1)
            elif len(set(df.phecode)) == 1:
                sign.append(i[:-1])
            else:
                sign.append(0)
        else:
            sign.append(1)
    return sign

icd9['icd'] = icd9.code.str.replace('.', '')
icd9['sign'] = getSign(icd9)
icd10['icd'] = icd10.code.str.replace('.', '')