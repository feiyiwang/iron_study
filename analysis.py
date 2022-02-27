import pandas as pd
import numpy as np
import tqdm

# All the necessary paths
# path to data in red_library
freq_path = 'finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9.afreq'
event_path = 'finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_detailed_longitudinal.txt.gz'
pca_path = 'finngen/library-red/finngen_R9/pca_1.0/data/finngen_R9.eigenvec.txt'
sex_path = 'finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_minimum_1.0.txt.gz'
# path to uploaded data
snp_path = 'finngen/green/FeiyiWang/all_variants_and_proxies.csv'
phecode_path = 'finngen/green/FeiyiWang/phecode_map.csv'


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
# shape of exposure_matrix is 392649 by 2460
# 392649 rows -> 392649 individuals collected in R9
# 2460 cols -> 1 finngen_id + 1 sex + 2458 dosages
exposure_matrix.to_csv('exposure_matrix.csv', index=None)
# delete dataframes to release some memory
del ped, freq, snp, demo


# Convert ICD codes to a list of phecodes
# load data
events = pd.read_csv(event_path, sep='\t')
phecode_map = pd.read_csv(phecode_path)

