import pandas as pd
import numpy as np
import tqdm

# paths
freq_path = 'finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9.afreq'
snp_path = 'finngen/green/FeiyiWang/all_variants_and_proxies.csv'
event_path = 'finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_detailed_longitudinal.txt.gz'
pca_path = 'finngen/library-red/finngen_R9/pca_1.0/data/finngen_R9.eigenvec.txt'
sex_path = 'finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_minimum_1.0.txt.gz'


# Obtain more accurate allele frequencies
freq = pd.read_csv(freq_path, sep='\t')
snp = pd.read_csv(snp_path)

snp = snp.merge(freq[['ID','ALT_FREQS']], left_on='sandbox_format', right_on='ID')
snp = snp.rename(columns={'ALT_FREQS':'sandbox_af'})
snp = snp.drop(columns=['Unnamed: 0', 'ID'])

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
ped = pd.read_csv('selected_snp.ped', sep=' ', header=None).rename(columns=ped_head)

# 3. build an exposure matrix
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
exposure_matrix.to_csv('exposure_matrix.csv', index=None)
# delete dataframes to release some memory
del ped, freq, snp, demo


# Convert ICDs to a list of phecodes
events = pd.read_csv(event_path, sep='\t')

