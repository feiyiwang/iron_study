import statsmodels.api as sm
from analysis_utils import *

# All the necessary paths
# path to data in red_library
freq_path = '/finngen/library-red/finngen_R9/genotype_plink_1.0/data/finngen_R9.afreq'
event_path = '/finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_detailed_longitudinal_1.0.txt.gz'
pca_path = '/finngen/library-red/finngen_R9/pca_1.0/data/finngen_R9.eigenvec.txt'
sex_path = '/finngen/library-red/finngen_R9/phenotype_1.0/data/finngen_R9_minimum_1.0.txt.gz'
# path to uploaded data
snp_path = '/finngen/green/FeiyiWang/all_variants_and_proxies.csv'
# source: R package - PheWAS
phecode_path = '/finngen/green/FeiyiWang/phecode_map.csv'
# source: https://phewascatalog.org/phecodes_icd10
phecode_path2 = '/finngen/green/FeiyiWang/phecode_icd10.csv'
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
del ped, freq, snp, demo


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
del pca


# Convert ICD codes to a list of phecodes
# 1. clean the event data
# load data
events = pd.read_csv(event_path, sep='\t')
# remove irrelevant codes
events = events[~events.SOURCE.isin(['OPER_OUT', 'OPER_IN'])]
events = events[(events.CATEGORY.str.startswith('ICD')) | (events.CATEGORY.str.match('\d+'))]
# keep only icd9 and icd10
events = events[events.ICDVER.isin(['9', '10'])][['FINNGENID', 'SOURCE', 'EVENT_AGE', 'CODE1']]
# remove n/a value
events = events[~events.CODE1.isna()]
# remove duplicates
events = events[~events.duplicated(subset=['FINNGENID', 'CODE1'])]
# split the data by icd version
events_icd9 = events[events.ICDVER == '9'][['FINNGENID', 'EVENT_AGE', 'CODE1']]
events_icd10 = events[events.ICDVER == '10'][['FINNGENID', 'EVENT_AGE', 'CODE1']]
# replace codes ending with lower case with upper case
to_fix = events_icd9[events_icd9.CODE1.str[-1].isin(['a', 'b', 'x'])].CODE1
for i, item in to_fix.iteritems():
    events_icd9.loc[i, 'CODE1'] = item.upper()

# 2. process the mapping data
# load data
phecode_map = pd.read_csv(phecode_path)
phecode_map2 = pd.read_csv(phecode_path2)  # icd 10 only
phecode_map2 = phecode_map2.rename(columns={'ICD10': 'code', 'PheCode': 'phecode'})
# 1) process the icd 9 mapping data
icd9 = phecode_map[phecode_map.vocabulary_id == 'ICD9CM'][['code', 'phecode']]
icd9['icd'] = icd9.code.str.replace('.', '')  # remove dot from the codes
icd9['sign'] = get_signs(icd9)
# 2) process the icd 10 mapping data
icd10 = phecode_map[phecode_map.vocabulary_id == 'ICD10CM'][['code', 'phecode']]
# some codes in icd10 contain '-'  e.g. A00-A09, T51-T65
icd10_ = icd10[icd10.code.str.contains('-')]  # len = 121
map_dict = {}
for i, row in icd10_.iterrows():
    code_type = row.code[0]
    code_range = [int(i[1:]) for i in row.code.split('-')]
    code_list = [str(i) for i in list(range(code_range[0], code_range[1] + 1))]
    code_list = ['0' + i if len(i) == 1 else i for i in code_list]
    for num in code_list:
        map_dict[code_type + num] = row.phecode
icd10_ = pd.DataFrame({'code': map_dict.keys(), 'phecode': map_dict.values()})
icd10 = pd.concat([icd10, icd10_, phecode_map2])
icd10 = icd10[(~icd10.phecode.isna()) & (~icd10.code.str.contains('-'))]
icd10 = icd10[~icd10.duplicated()]
icd10['icd'] = icd10.code.str.replace('.', '')
icd10['sign'] = get_signs(icd10)
# split to two parts by checking duplicates - one phecode to many icd codes
dup_list = icd10[icd10.icd.duplicated()].icd.tolist()
icd10_dup = icd10[icd10.code.isin(dup_list)]  # len = 8110  7457
icd10 = icd10[~icd10.code.isin(dup_list)]  # len = 66722  66909
icd10_dup = icd10_dup.sort_values(by='icd')
icd10 = icd10.sort_values(by='icd')

# 3. ICD 9 mapping
# 1) exact mapping
events_icd9 = events_icd9.merge(icd9[['icd', 'phecode']], 'left', left_on='CODE1', right_on='icd')
print_result(events_icd9)  # 0.0%
# 2) mapping by the first 4 elements
# if all the codes with the same first 4 elements have the same phecode
icd9_ = icd9[~icd9.sign.isin([0, 1])][['phecode', 'sign']]
icd9_ = icd9[~icd9_.duplicated()]
events_icd9_ = get_events_subset(events_icd9)
events_icd9_ = events_icd9_.merge(icd9[['sign', 'phecode']], 'left', left_on='code_', right_on='sign')
events_icd9_ = events_icd9_[~events_icd9_.phecode.isna()]
events_icd9.loc[events_icd9_['index'].tolist(), 'phecode'] = events_icd9_.phecode.tolist()
print_result(events_icd9)  # 29.2%
# if these 5-digit Finnish codes can be converted to standardized 4-digit codes only
events_icd9_ = get_events_subset(events_icd9)
events_icd9_['phecode'] = get_phecodes(events_icd9_, icd9)
events_icd9_ = events_icd9_[~events_icd9_.phecode.isin(['fail', 'na'])]
events_icd9.loc[events_icd9_['index'].tolist(), 'phecode'] = events_icd9_.phecode.tolist()
print_result(events_icd9)  # 68.7%

# 4. ICD 10 mapping
events_icd10_dup = events_icd10[events_icd10.CODE1.isin(dup_list)]
events_icd10 = events_icd10[~events_icd10.CODE1.isin(dup_list)]
''' For events_icd10 '''
# 1) exact mapping
# if each icd code is mapped to only one phecode
events_icd10 = events_icd10.merge(icd10[['icd', 'phecode']], 'left', left_on='CODE1', right_on='icd')
print_result(events_icd10)  # 77.2%
# if the code if from those A00-A09, T51-T65, etc.
events_icd10_ = events_icd10[events_icd10.phecode.isna()]
events_icd10_['code_'] = events_icd10_.CODE1.str[:3]
events_icd10_ = events_icd10_[['FINNGENID', 'code_']]
events_icd10_['index'] = events_icd10_.index
events_icd10_ = events_icd10_.merge(icd10_, 'left', left_on='code_', right_on='code')
events_icd10_ = events_icd10_[~events_icd10_.phecode.isna()]
events_icd10.loc[events_icd10_['index'].tolist(), 'phecode'] = events_icd10_.phecode.tolist()
print_result(events_icd10)  # 82.3%
# 2) mapping by the first n-1 elements in the code
# e.g. Given O9980 in Finngen, search O998 and check if all codes begin with O998 are mapped to the same phecode
events_icd10_ = events_icd10[events_icd10.phecode.isna()][['FINNGENID', 'CODE1']]
events_icd10_['index'] = events_icd10_.index
code = events_icd10_.CODE1.unique().tolist()
phecode = [get_phecodes(i[:-1], icd10) for i in tqdm.tqdm(code)]
events_icd10_ = events_icd10_.merge(pd.DataFrame({'CODE1':code, 'phecode':phecode}), 'left')
events_icd10_ = events_icd10_[~events_icd10_.phecode.isin('?')]
events_icd10.loc[events_icd10_['index'].tolist(), 'phecode'] = events_icd10_.phecode.tolist()
print_result(events_icd10)  # 91.0%
events_icd10 = events_icd10[['FINNGENID', 'CODE1', 'phecode']]
''' For events_icd10_dup '''
# exact mapping
# merge a new dataset to concat the events_icd10
events_icd10_ = events_icd10_dup[~events_icd10_dup.duplicated()][['FINNGENID', 'CODE1']]
events_icd10_ = events_icd10_.merge(icd10_dup[['icd', 'phecode']], left_on='CODE1', right_on='icd')
events_icd10 = pd.concat([events_icd10, events_icd10_], axis=0)
print_result(events_icd10)  # 92.9%

# 5. merge events and clean the data
# icd10 accounts for 95.95% of all the codes in events
events = pd.concat([events_icd10, events_icd9], axis=0)
print_result(events) # 92.0%
events.to_csv('events.csv', index=None)
# get a list of unique phecodes for modeling
phecode_list = sorted(phecode_map.phecode.unique().tolist())
events = events[events.phecode.isin(phecode_list)][['FINNGENID', 'phecode']]
events = events.rename(columns={'FINNGENID':'finngen_id'})


# Convert events to a matrix of outcomes
confounding = pd.read_csv('confounding.csv')
outcome_matrix = confounding[['finngen_id']]
n_cases = []
# add phecodes to the dataframe
for i in tqdm.tqdm(range(len(phecode_list))):
    id_list = events[events.phecode == phecode_list[i]].finngen_id.tolist()
    case_bool = np.select([(outcome_matrix.finngen_id.isin(id_list)),
                         (~outcome_matrix.finngen_id.isin(id_list))], [1, 0])
    outcome_matrix[phecode_list[i]] = case_bool
    n_cases.append(len(id_list))
# shape of outcome_matrix is 392649 by 1861
# 392649 rows -> 392649 individuals collected in R9
# 2459 cols -> 1 finngen_id + 1860 PheCodes
outcome_matrix.to_csv('outcome_matrix.csv', index=None)
# delete dataframes to release some memory
del events


# Build a for-loop for logistic regression
# 1. prepare a phecode definition table for sex and n_case
phecode_path3 = '/finngen/green/FeiyiWang/phecode_sex.csv'
phecode_sex = pd.read_csv(phecode_path3) # 166 sex-specified
phecode_def = pd.DataFrame({'phecode':phecode_list, 'n_cases':n_cases})
phecode_def = phecode_def[phecode_def.n_cases >= 100]
phecode_def = phecode_def.merge(phecode_sex, 'left') # 147 sex-specified
phecode_def['sex'] = np.select([(phecode_def.sex == 'Female'),
                                (phecode_def.sex == 'Male'),
                                (phecode_def.sex.isna())
                                ], [1, 0, -1])

# shape of phecode_def is 1478 by 3
# 1478 rows -> 1478 Phecodes
# 3 cols -> 1 phecode_id + 1 n_cases + 1 sex
phecode_def.to_csv('phecode_def.csv', index=None)

# 2. clean all the matrices for modeling
# remove 143 individuals without sex & age
ids_to_remove = confounding[confounding.sex == ''].finngen_id.tolist()
if ids_to_remove == confounding[confounding.age.isna()].finngen_id.tolist():
    confounding = confounding[~confounding.finngen_id.isin(ids_to_remove)]
    exposure_matrix = exposure_matrix[~exposure_matrix.finngen_id.isin(ids_to_remove)]
    outcome_matrix = outcome_matrix[~outcome_matrix.finngen_id.isin(ids_to_remove)]

# fix 15146 individuals without pca information. Replace them with 0s
confounding = confounding.fillna(0.0)
# check the distribution of confounding

ax = confounding.plot.hist(column=['age'], alpha=0.5, figsize=(10, 8))

stats = pd.DataFrame(columns=['Coef.','Std.Err.','z','P>|z|','exposure','outcome','n_cases', 'n_cohort'])
for snp in tqdm.tqdm(exposure_matrix.columns[1:]):
    x = pd.concat([exposure_matrix[[snp]], confounding.iloc[:, 1:]], axis=1)
    x = sm.add_constant(x)
    for _, row in phecode_def.iterrows():
        y = outcome_matrix[[row.phecode]]
        if type(row.sex) == int:
            x = x[x.sex == row.sex]
            x = x.drop(columns=['sex'])
            y = y[y.index.isin(x.index)]
        model = sm.Logit(y, x).fit()
        stat = model.summary2().tables[1].loc[snp, ['Coef.', 'Std.Err.', 'z', 'P>|z|']]
        stat['exposure'] = snp
        stat['outcome'] = row.phecode
        stat['n_cases'] = row.n_cases
        stat['n_cohort'] = len(y)
        stats = stats.append(stat)
