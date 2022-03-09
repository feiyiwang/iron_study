import tqdm
import pandas as pd
import numpy as np
from scipy import stats


def get_signs(data):
    sign = []
    for i in tqdm.tqdm(data.icd):
        if len(i) == 5:
            df = data[data.icd.str.startswith(i[:4])]
            if len(df) == 1:
                sign.append(1)
            elif len(set(df.phecode)) == 1:
                sign.append(i[:4])
            else:
                sign.append(0)
        else:
            sign.append(1)
    return sign


def get_events_subset(data):
    data_ = data[(data.phecode.isna()) & (data.CODE1.str.len() == 5)]
    data_['code_'] = data_.CODE1.str[:4]
    data_ = data_[['FINNGENID', 'code_']]
    data_['index'] = data_.index  # index will change after merging
    return data_


def get_phecodes(data, icd_map):
    code = []
    for i in tqdm.tqdm(data.code_):
        try:
            next_code = icd_map.loc[icd_map[icd_map.icd == i].index + 1, 'icd'].values[0]
            if (len(next_code) != 5) or (next_code[:4] != i):
                this_phecode = icd_map[icd_map.icd == i].phecode.values[0]
                code.append(this_phecode)
            else:  # this code need further discussion
                code.append('fail')
        except IndexError:  # this code_ does not exist in phecode_map
            code.append('na')
    return code


def get_phecode(icd, icd_map):
    sub_map = icd_map[icd_map.icd.str.startswith(icd)]
    if len(sub_map) == 0:
        return '?'
    elif len(set(sub_map.phecode)) == 1:
        return sub_map.phecode.tolist()[0]
    else:
        return '?'


def print_result(data):
    print(
        len(data[~data.phecode.isna()]), '/',
        len(data), '=',
        len(data[~data.phecode.isna()]) / len(data)
    )


def get_summary(model, X, y, var_names):
    """
    :param model: a well-trained model instance
    :param X: a DataFrame of inputs
    :param y: a Series of output
    :param var_names: a list of confounding names
    :return: a DataFrame of key statistics
    """
    coef = np.append(model.intercept_, model.coef_)
    X_matrix = np.append(np.ones((len(X), 1)), X, axis=1)
    y_hat = model.predict(X)
    degree_of_freedom = X_matrix.shape[0] - X_matrix.shape[1]
    MSE = (sum((y_hat - y) ** 2)) / degree_of_freedom
    var_coef = MSE * (np.linalg.inv(np.dot(X_matrix.T, X_matrix)).diagonal())
    std_err_coef = np.sqrt(var_coef)
    t_stat_coef = coef / std_err_coef
    p_values_coef = [2 * (1 - stats.t.cdf(np.abs(i), degree_of_freedom)) for i in t_stat_coef]
    t_half_alpha = stats.t.ppf(1 - 0.025, degree_of_freedom)
    ci1_coef = [beta - t_half_alpha * se_beta for beta, se_beta in zip(coef, std_err_coef)]
    ci2_coef = [beta + t_half_alpha * se_beta for beta, se_beta in zip(coef, std_err_coef)]
    summary_df = pd.DataFrame({
        'coef': np.round(coef, 4),
        'std_err': np.round(std_err_coef, 3),
        't_stat': np.round(t_stat_coef, 3),
        'p_value': np.round(p_values_coef, 3),
        'ci_1': np.round(ci1_coef, 3),
        'ci_2': np.round(ci2_coef, 3),
    }, index=['const'] + var_names)
    return summary_df

