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
phecode_sex = pd.read_csv('/finngen/green/FeiyiWang/phecode_sex.csv') # 166 sex-specified


# 2. Remove individuals without PCA info
# finngen ids are in the same order across the datasets
confounding = confounding[confounding.PC1 != 0]
confounding = confounding.drop(columns=['finngen_id'])
exposure_matrix = exposure_matrix[snp_final.sandbox_format.to_list()]
exposure_matrix = exposure_matrix[exposure_matrix.index.isin(confounding.index)]
outcome_matrix = outcome_matrix[outcome_matrix.index.isin(confounding.index)]
outcome_matrix = outcome_matrix.drop(columns=['id'])
# replace N/A values with False in outcome_matrix
outcome_matrix = outcome_matrix.fillna(False)


# 3. Prepare a phecode definition table for sex and n_case
phecode_def = pd.DataFrame({
    'phecode': outcome_matrix.columns,
    'n_cases': outcome_matrix.sum(axis=0)
})
phecode_def['phecode'] = phecode_def.phecode.astype(float)
phecode_def['phe_id'] = phecode_def.index
phecode_def = phecode_def.merge(phecode_sex, 'left')
phecode_def['sex'] = np.select([(phecode_def.sex == 'Female'),
                                (phecode_def.sex == 'Male'),
                                (phecode_def.sex.isna())
                                ], [1, 0, -1])  # 1836 rows x 4 columns
# remove phecodes that has fewer than 100 individuals
phecode_def = phecode_def[phecode_def.n_cases >= 100] # 147 sex-specified  # 1215 rows x 4 columns


# 4. Build a for-loop for logistic regression
start = datetime.datetime.now()
results_df = pd.DataFrame(columns=['Coef.', 'Std.Err.', 'z', 'P>|z|', 'outcome', 'n_cases', 'n_cohort'])
for snp in tqdm.tqdm(exposure_matrix.columns):
    inputs = pd.concat([exposure_matrix[[snp]], confounding], axis=1)
    inputs = sm.add_constant(inputs)

    def modeling(subset, x=inputs):
        # subset column order:
        # 0-phecode-float,
        # 1-n_cases-int,
        # 2-phe_id-str,
        # 3-sex-int
        y = outcome_matrix[[subset[2]]]
        if subset[3] > -1:
            x = x[x.sex == subset[3]]
            x = x.drop(columns=['sex'])
            y = y[y.index.isin(x.index)]
        try:
            model = sm.Logit(y, x).fit()
            stat = model.summary2().tables[1].loc[snp, ['Coef.', 'Std.Err.', 'z', 'P>|z|']]
        except PerfectSeparationError:
            # When event per variable is too low,
            # the logistic model is under the risk of complete separation or quasi⁃complete separation.
            # https://m.medsci.cn/article/show_article.do?id=606319405019
            stat = pd.Series({'Coef.': 0., 'Std.Err.': 0., 'z': 0., 'P>|z|': 0.}, name=snp)
        except np.linalg.LinAlgError:
            stat = pd.Series({'Coef.': 1., 'Std.Err.': 1., 'z': 1., 'P>|z|': 1.}, name=snp)
        stat['outcome'] = subset[2]
        stat['n_cases'] = subset[1]
        stat['n_cohort'] = len(y)
        return stat

    pool = mp.Pool(processes=(mp.cpu_count() - 1)) # mp.cpu_count() = 16 in sandbox
    results = pool.map(modeling, phecode_def.to_numpy())
    pool.close()
    pool.join()
    for res in tqdm.tqdm(results):
        results_df = results_df.append(res)
end = datetime.datetime.now()
print(end - start)
# ~52 mins for one snp