import datetime
import multiprocessing as mp
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
# add cols to events according to createPhenotypes requirements
events['vocabulary_id'] = np.select([(events.ICDVER == '9'), (events.ICDVER == '10')], ['ICD9CM', 'ICD10CM'])
events = events.rename(columns={'FINNGENID': 'id', 'CODE1': 'code'})
events['count1'] = 1
events = events[['id', 'vocabulary_id', 'code', 'count1']]
events.to_csv('events.csv', index=None)
# create a sex_df for createPhenotype() in R
sex_df = confounding[confounding.sex != ''][['finngen_id', 'sex']]
sex_df['sex'] = np.select([(sex_df.sex == 1), (sex_df.sex == 0)], ['M', 'F'])
sex_df = sex_df.rename(columns={'finngen_id': 'id'})
sex_df.to_csv('sex_df.csv', index=None)
# the rest will be done in R


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
# remove phecodes that has fewer than 100 individuals
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

# remove 15146 individuals without pca information
confounding = confounding[~confounding.isna()]#.fillna(0.0)
# check the distribution of confounding

ax = confounding.plot.hist(column=['age'], alpha=0.5, figsize=(10, 8))

stats = pd.DataFrame(columns=['Coef.','Std.Err.','z','P>|z|','outcome','n_cases', 'n_cohort'])
for snp in tqdm.tqdm(exposure_matrix.columns[1:]):
    x = pd.concat([exposure_matrix[[snp]], confounding.iloc[:, 1:]], axis=1)
    x = sm.add_constant(x)
    for _, row in phecode_def.iterrows():
        y = outcome_matrix[[row.phecode]]
        if type(row.sex) > -1:
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
# 1:26:58 for one snp


start = datetime.datetime.now()
results_df = pd.DataFrame()
for snp in tqdm.tqdm(exposure_matrix.columns[67:68]):
    inputs = pd.concat([exposure_matrix[[snp]], confounding.iloc[:, 1:]], axis=1)
    inputs = sm.add_constant(inputs)

    def modeling(subset, x=inputs):
        y = outcome_matrix[[subset[0]]]
        if subset[2] > -1:
            x = x[x.sex == subset[2]]
            x = x.drop(columns=['sex'])
            y = y[y.index.isin(x.index)]
        model = sm.Logit(y, x).fit()
        stat = model.summary2().tables[1].loc[snp, ['Coef.', 'Std.Err.', 'z', 'P>|z|']]
        stat['outcome'] = subset[0]
        stat['n_cases'] = subset[1]
        stat['n_cohort'] = len(y)
        return stat

    pool = mp.Pool(processes=(mp.cpu_count() - 1))
    results = pool.map(modeling, phecode_def.to_numpy())
    pool.close()
    pool.join()
    results_df = pd.concat(results)
end = datetime.datetime.now()
print(end - start)
# 1:02:29 for one snp


from scipy.stats import uniform
from scipy.stats import randint
import numpy as np
import matplotlib.pyplot as plt

# sample data
df = pd.DataFrame({'gene' : ['gene-%i' % i for i in np.arange(10000)],
'pvalue' : uniform.rvs(size=10000),
'chromosome' : ['ch-%i' % i for i in randint.rvs(0,12,size=10000)]})

# -log_10(pvalue)
df['minuslog10pvalue'] = -np.log10(df.pvalue)
df.chromosome = df.chromosome.astype('category')
df.chromosome = df.chromosome.cat.set_categories(['ch-%i' % i for i in range(12)], ordered=True)
df = df.sort_values('chromosome')

# How to plot gene vs. -log10(pvalue) and colour it by chromosome?
df['ind'] = range(len(df))
df_grouped = df.groupby(('chromosome'))

# manhattan plot
fig = plt.figure(figsize=(14, 8)) # Set the figure size
ax = fig.add_subplot(111)
colors = ['darkred','darkgreen','darkblue', 'gold']
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=colors[num % len(colors)], ax=ax)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)

# set axis limits
ax.set_xlim([0, len(df)])
ax.set_ylim([0, 3.5])

# x axis label
ax.set_xlabel('Chromosome')

# show the graph
plt.show()

# fam id
# robust error
#
# gene name vs rsid




