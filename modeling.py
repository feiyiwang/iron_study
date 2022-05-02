import tqdm
import pandas as pd
import numpy as np
import datetime
import multiprocessing as mp
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError

# 1. Load all the datasets
exposure_matrix = pd.read_csv('exposure_matrix.csv')
confounding = pd.read_csv('confounding.csv')
outcome_matrix = pd.read_csv('outcome_matrix.csv')
snp_final = pd.read_csv('/finngen/green/FeiyiWang/all_variants_and_proxies_finn_af_final.csv')


# 2. Remove individuals without PCA info
# finngen ids are in the same order across the datasets
confounding = confounding[confounding.PC1 != 0]
exposure_matrix = exposure_matrix[snp_final.sandbox_format.to_list()]
exposure_matrix = exposure_matrix[exposure_matrix.index.isin(confounding.index)]
outcome_matrix = outcome_matrix[outcome_matrix.id.isin(confounding.finngen_id)]
outcome_matrix = outcome_matrix.rename(columns={'id': 'finngen_id'})
# # replace N/A values with False in outcome_matrix
# outcome_matrix = outcome_matrix.fillna(False)


# 3. Prepare a phecode definition table for sex and n_case
phecode_def = (len(outcome_matrix) - outcome_matrix.isna().sum()).to_frame('n_cohort')
phecode_def['phecode'] = phecode_def.index,
phecode_def['n_cases'] = outcome_matrix.sum(axis=0)
phecode_def = phecode_def.iloc[1:, :]
# remove phecodes that has fewer than 100 individuals
phecode_def = phecode_def[phecode_def.n_cases >= 100]
# calculate n_males and n_females
m, f = [], []
for i in tqdm.tqdm(phecode_def.phecode):
    df = outcome_matrix[['finngen_id', i]].merge(confounding[['finngen_id', 'sex']], 'left', on='finngen_id')
    sex_proportion = df[df[i] == True].sex.value_counts()
    n_males = sex_proportion[0] if 0 in sex_proportion.keys() else 0
    n_females = sex_proportion[1] if 1 in sex_proportion.keys() else 0
    m.append(n_males)
    f.append(n_females)

phecode_def['n_males'] = m
phecode_def['n_females'] = f
phecode_def['sex'] = np.select([
    (phecode_def.n_males == 0), (phecode_def.n_females == 0),
    ((phecode_def.n_males != 0) | (phecode_def.n_females != 0))
], [1, 0, -1])
# 33 female-specified  # 7 male-specified


# 4. Clean the datasets accordingly
exposure_matrix['finngen_id'] = confounding.finngen_id
exposure_matrix = exposure_matrix.sort_values(by='finngen_id')
confounding = confounding.sort_values(by='finngen_id')
outcome_matrix = outcome_matrix[phecode_def.phecode.tolist()+['finngen_id']].sort_values(by='finngen_id')
assert len(exposure_matrix) == len(outcome_matrix) == len(confounding)
exposure_matrix.index = range(len(exposure_matrix))
confounding.index = range(len(exposure_matrix))
outcome_matrix.index = range(len(exposure_matrix))


# 5. Build a for-loop for logistic regression
start = datetime.datetime.now()
results_df = pd.DataFrame(columns=['Coef.', 'Std.Err.', 'z', 'P>|z|', 'outcome', 'n_cases', 'n_cohort'])
for snp in tqdm.tqdm(exposure_matrix.columns[:-1]):
    inputs = pd.concat([exposure_matrix[[snp]], confounding.iloc[:, 1:]], axis=1)
    inputs = sm.add_constant(inputs)

    def modeling(subset, x=inputs):
        # subset column order:
        # 0-n_cohort-int,
        # 1-n_cases-int,
        # 2-phe_id-str,
        # 3-sex-int
        y = outcome_matrix[[subset[2]]][~outcome_matrix[subset[2]].isna()]
        x = x[x.index.isin(y.index)]
        if subset[3] > -1:
            x = x[x.sex == subset[3]]
            x = x.drop(columns=['sex'])
            y = y[y.index.isin(x.index)]
        try:
            model = sm.Logit(y, x).fit(disp=0)
            stat = model.summary2().tables[1].loc[snp, ['Coef.', 'Std.Err.', 'z', 'P>|z|']]
        except PerfectSeparationError:
            # When event per variable is too low,
            # the logistic model is under the risk of complete separation or quasi⁃complete separation.
            # https://m.medsci.cn/article/show_article.do?id=606319405019
            num = 2.
            stat = pd.Series({'Coef.': num, 'Std.Err.': num, 'z': num, 'P>|z|': num}, name=snp)
            print('PerfectSeparationError: ', subset[2])
        except np.linalg.LinAlgError:
            num = 3.
            stat = pd.Series({'Coef.': num, 'Std.Err.': num, 'z': num, 'P>|z|': num}, name=snp)
            print('LinAlgError: ', subset[2])
        except ValueError:
            # ValueError: Pandas data cast to numpy dtype of object. Check input data with np.asarray(data).
            num = 4.
            stat = pd.Series({'Coef.': num, 'Std.Err.': num, 'z': num, 'P>|z|': num}, name=snp)
            print('ValueError: ', subset[2])
        except Exception as e:
            num = 5.
            stat = pd.Series({'Coef.': num, 'Std.Err.': num, 'z': num, 'P>|z|': num}, name=snp)
            print(e, ': ', subset[2])
        stat['outcome'] = subset[2]
        stat['n_cases'] = subset[1]
        stat['n_cohort'] = subset[0]
        return stat

    pool = mp.Pool(processes=(mp.cpu_count() - 1)) # mp.cpu_count() = 16 in sandbox
    results = pool.map(modeling, phecode_def.to_numpy())
    pool.close()
    pool.join()
    for res in results:
        results_df = results_df.append(res)
end = datetime.datetime.now()
print(end - start)
# ~17 mins for one snp
# ~39 hours for the whole loop