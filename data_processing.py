from utils import *

# All the necessary paths
# path to data in red_library
freq_path = '/finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9.afreq'
event_path = '/finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_detailed_longitudinal_1.0.txt.gz'
pca_path = '/finngen/library-red/finngen_R9/pca_1.0/data/finngen_R9.eigenvec.txt'
sex_path = '/finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_minimum_1.0.txt.gz'
# path to the uploaded data
snp_path = '/finngen/green/FeiyiWang/all_variants_and_proxies.csv'
# source: R package - PheWAS
phecode_path = '/finngen/green/FeiyiWang/phecode_map.csv'
# source: https://phewascatalog.org/files/phecode_definitions1.2.csv.zip
phecode_path3 = '/finngen/green/FeiyiWang/phecode_sex.csv'


# Obtain more accurate allele frequencies
# load data
freq = pd.read_csv(freq_path, sep='\t')
snp = pd.read_csv(snp_path)
# left join snp and freq
snp = snp.merge(freq[['ID', 'ALT_FREQS']], left_on='sandbox_format', right_on='ID')
snp = snp.rename(columns={'ALT_FREQS': 'sandbox_af'})
# drop unnecessary cols
snp = snp.drop(columns=['Unnamed: 0', 'ID'])
# save the updated snp dataframe
snp.to_csv('all_variants_and_proxies.csv', index=None)


# Build a matrix of genetic dosages
# 1. save a list of SNP positions for obtaining genetic data at individual level via plink
snp_sub = snp[(~snp.sandbox_format.isna()) & (snp.match == True) & (~snp.sandbox_format.duplicated())]
snp_list = snp_sub.sandbox_format.tolist()
snp_strings = ''
for i in snp_list:
    snp_strings += i + ','

# 2. run plink to obatin a .ped file directly in notebook
# !plink --bfile /finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9 --snps $snp_strings --recode --out selected_snp
ped_head = {
    0: 'family_id',
    1: 'finngen_id',
    2: 'father_id',
    3: 'mother_id',
    4: 'sex',
    5: 'phenotype'
}
# In .ped file, the first 6 cols are listed as ped_head
# and then each two cols in the rest are a pair of alleles tied to a specific SNP
ped = pd.read_csv('selected_snp.ped', sep=' ', header=None).rename(columns=ped_head)

# 3. build an exposure matrix
# get a list of reference alleles
snp_ref = snp_sub.finn_ref.tolist()
# set up a dataframe with ids for the dosage matrix
exposure_matrix = ped[['finngen_id']]
# add dosages to the dataframe
for i in tqdm.tqdm(range(len(snp_list))):
    allele1 = np.select([(ped[6 + 2 * i] == snp_ref[i]), (ped[6 + 2 * i] != snp_ref[i])], [1, 0])
    allele2 = np.select([(ped[6 + 2 * i + 1] == snp_ref[i]), (ped[6 + 2 * i + 1] != snp_ref[i])], [1, 0])
    dosage = allele1 + allele2
    exposure_matrix[snp_list[i]] = dosage
# shape of exposure_matrix is 392649 by 2459
# 392649 rows -> 392649 individuals collected in R9
# 2459 cols -> 1 finngen_id + 2458 SNPs
exposure_matrix.to_csv('exposure_matrix.csv', index=None)
# delete dataframes to release some memory
del ped, freq, snp


# Build a matrix of confounding
confounding = exposure_matrix[['finngen_id']]
# add sex to the dataframe
demo = pd.read_csv(sex_path, sep='\t')
demo = demo.rename(columns={'FINNGENID': 'finngen_id', 'SEX': 'sex', 'BL_AGE': 'age'})
confounding = confounding.merge(demo[['finngen_id', 'sex']], 'left')
# 1 is female; 0 is male
confounding['sex'] = np.select([(confounding.sex == 'female'),
                                (confounding.sex == 'male'),
                                (confounding.sex.isna())], [1, 0, ''])
# add 10 principal components of ancestry to the data
pca = pd.read(pca_path)
confounding = confounding.merge(pca.iloc[:, 1:12], 'left', left_on='finngen_id', right_on='IID')
confounding = confounding.drop(columns=['IID'])
# add age at the measurement to the data
confounding = confounding.merge(demo[['finngen_id', 'age']], 'left')
# shape of exposure_matrix is 392649 by 13
# 392649 rows -> 392649 individuals collected in R9
# 13 cols -> 1 finngen_id + 1 sex + 10 principal components + 1 age
confounding.to_csv('confounding.csv', index=None)
# delete dataframes to release some memory
del pca, demo


# Convert ICD codes to a list of phecodes
# load data
events = pd.read_csv(event_path, sep='\t')
# remove irrelevant codes
events = events[(events.CATEGORY.str.startswith('ICD')) | (events.CATEGORY.str.match('\d+'))]
# keep only icd9 and icd10
events = events[events.ICDVER.isin(['9', '10'])][['FINNGENID', 'SOURCE', 'CODE1', 'ICDVER']]
# remove n/a value
events = events[~events.CODE1.isna()]
# add dot back to ICD codes
events['code_len'] = events.CODE1.str.len()
events_1 = events[events.code_len == 3]
events_1['code'] = events_1.CODE1
events_2 = events[events.code_len > 3]
events_2['code'] = events_2.CODE1.str[:3]+'.'+events_2.CODE1.str[3:]
events = pd.concat([events_1, events_2], axis=0)
# add cols to events according to createPhenotypes requirements
events['vocabulary_id'] = np.select([(events.ICDVER == '9'), (events.ICDVER == '10')], ['ICD9CM', 'ICD10CM'])
events = events.rename(columns={'FINNGENID': 'id'})
events['count1'] = 1
events = events[['id', 'vocabulary_id', 'code', 'count1']]
events.to_csv('events.csv', index=None)
# create a sex_df for createPhenotype() in R
sex_df = confounding[confounding.sex != ''][['finngen_id', 'sex']]
sex_df['sex'] = np.select([(sex_df.sex == 1), (sex_df.sex == 0)], ['F', 'M'])
sex_df = sex_df.rename(columns={'finngen_id': 'id'})
sex_df.to_csv('sex_df.csv', index=None)
# the rest will be done in R