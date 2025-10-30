#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as smm
import argparse
from scipy import stats
import warnings
import time
from multiprocessing import Pool, Process


# function to parse output of Model2 predictions
# Define a function to parse M2 output
def parse_M2output(M2_file_path, RepName, Condition):
    summary_data = pd.read_csv(M2_file_path, sep='\t', dtype={'contig': str, 'position': int, 'site': str})
    summary_data.rename(columns={'coverage': f'coverage_{Condition}', 'stoichiometry': f'stoichiometry_{Condition}',
                                 'probability': f'probability_{Condition}'}, inplace=True)
    return summary_data


def calculate_nested_f_statistic(small_model, big_model):
    # Given two fitted GLMs, the larger of which contains the parameter space of the smaller, return the F Stat and P value corresponding to the larger model adding explanatory power#
    addtl_params = big_model.df_model - small_model.df_model
    f_stat = (small_model.deviance - big_model.deviance) / (addtl_params * big_model.scale)
    df_numerator = addtl_params
    # use fitted values to obtain n_obs from model object:
    df_denom = (big_model.fittedvalues.shape[0] - big_model.df_model)
    p_value = stats.f.sf(f_stat, df_numerator, df_denom)
    return (f_stat, p_value)


def combine_data(df):
    # Create a dictionary to hold combined values for each measurement
    combined_data = {}
    print("\n\ncombining here\n\n",df)
    # Get unique measurements (coverage_KOPN, stoichiometry_KOPN, probability_KOPN, etc.)
    measurements = set('_'.join(col.split('_')[:2]) for col in df.columns if '_' in col)

    for measurement in measurements:
        # Select columns that match the current measurement
        measurement_cols = [col for col in df.columns if col.startswith(measurement)]

        # Combine values from matching columns for each row
        combined_data[measurement] = df[measurement_cols].apply(list, axis=1)

    # Create a new DataFrame from the combined data
    new_df = pd.DataFrame(combined_data)
    column_order = ['contig', 'position', 'site'] + sorted(list(measurements))
    # Add the non-measurement columns ('contig', 'position', 'site')
    new_df[['contig', 'position', 'site']] = df[['contig', 'position', 'site']]
    new_df = new_df[column_order]
    return new_df


def custom_round(num):
    if num % 1 >= 0.5:  # Check if the decimal part is 0.5 or greater
        num += 0.5  # If true, add 0.5 to the number
    return round(num)


def make_bernoulli_data(maindf):
    # Given a data frame containing stoichiometry and coverage each replicate of the condition to compare, recapitulate read_level outcome#
    temp_success = []
    coverage_columns = [col for col in maindf.index if col.startswith('coverage_')]
    # Create lists based on column values with varying lengths
    temp_condition = []
    for col in coverage_columns:
        prefix = col.split('_')[1]  # Extracting the list from the Series
        temp_condition.extend([prefix] * int(maindf[col].sum()))
        for success_val, failure_val in zip(maindf['success' + prefix], maindf['failure' + prefix]):
            temp_success.extend([1] * round(success_val))  # Append 1s based on 'successKOPN' values
            temp_success.extend([0] * round(failure_val))  # Extend the result list based on the length of each list
    data = {'Success': temp_success, 'Condition': temp_condition}
    return (pd.DataFrame(data))


def calculate_statistic(small_model, big_model):
    from scipy.stats import chi2
    """Like calculate_nested_f_statistic but corrects for overdispersion"""
    # deviance_diff = small_model.deviance - big_model.deviance
    # df_diff = big_model.df_model - small_model.df_model
    deviance_diff = small_model.deviance - big_model.deviance
    df_diff = big_model.df_model - small_model.df_model
    p_value = chi2.sf(deviance_diff, df_diff)
    return p_value


def test_group_effect_bernoulli(df):
    big_model = smf.glm('stoichiometry ~ condition', data=df,
                        family=sm.families.Binomial(), freq_weights=df['coverage']).fit()
    small_model = smf.glm('stoichiometry ~ 1', data=df,
                          family=sm.families.Binomial(), freq_weights=df['coverage']).fit()
    return calculate_statistic(small_model, big_model)


def test_group_effect_binom(df, conditions):
    select_columns = df.filter(regex='^success')

    # Concatenate values into a single list
    successes = [val for sublist in select_columns.values.flatten() for val in sublist]

    select_columns = df.filter(regex='^failure')

    # Concatenate values into a single list
    failures = [val for sublist in select_columns.values.flatten() for val in sublist]
    #print(len(conditions), len(successes), len(failures))
    # Creating a DataFrame with extracted values
    data = {'Condition': conditions, 'Success': successes, 'Failures': failures}
    df_data = pd.DataFrame(data)
    # print(df_data)
    df_data = df_data[~((df_data.Success == 0) & (df_data.Failures == 0))]
    # Fitting the GLM models
    big_model = smf.glm('Success + Failures ~ Condition', family=sm.families.Binomial(), data=df_data).fit()
    small_model = smf.glm('Success + Failures ~ 1', family=sm.families.Binomial(), data=df_data).fit()

    # Calculating F-statistic and p-value for nested models
    fstat, pvalue = calculate_nested_f_statistic(small_model, big_model)  # You need to define this function
    return fstat, pvalue


# Example usage:
# new_df = combine_data(your_original_df)

# function to check if all values in a list are NaN
def all_nan(lst):
    return all(np.isnan(x) for x in lst)


def calculate_average_stoichiometry(row, cond):
    coverage_col = f'coverage_{cond}'
    stoichiometry_col = f'stoichiometry_{cond}'

    coverage = np.array(row[coverage_col])
    stoichiometry = np.array(row[stoichiometry_col])

    # Calculate the sum of coverage
    sum_coverage = np.nansum(coverage)

    # Calculate the sum of element-wise multiplication of coverage and stoichiometry
    sum_elements = np.nansum(coverage * stoichiometry)

    if sum_coverage == 0:
        return np.nan
    else:
        return sum_elements / sum_coverage


# Function to check if only one value is non-zero in a list
def one_non_zero(arr):
    non_zero_count = sum(1 for i in arr if i != 0)
    return non_zero_count == 1


def combine_samples(args):
    df = pd.DataFrame()
    data_to_compare = pd.read_csv(args.data_file, sep='\t')
    for ind in data_to_compare.index:
        temp_df = parse_M2output(data_to_compare['M2_file_path'][ind], data_to_compare['RepName'][ind],
                                 data_to_compare['Condition'][ind])
        if df.empty:
            df = temp_df
        else:
            df = df.merge(temp_df, on=['contig', 'position', 'site'], how='outer')
    df = combine_data(df)
    for cond in data_to_compare['Condition'].unique():
        df = df[~(df["coverage_" + cond].apply(all_nan)) & (
            df['probability_' + cond].apply(lambda x: any(val > args.M2_threshold for val in x))) \
                & (df['stoichiometry_' + cond].apply(lambda x: any(val > args.Stoichiometry_threshold for val in x)))]
        df['coverage_' + cond] = df['coverage_' + cond].apply(lambda x: np.nan_to_num(x))
        df['stoichiometry_' + cond] = df['stoichiometry_' + cond].apply(lambda x: np.nan_to_num(x))
        df[f'avg_stoichi_{cond}'] = df.apply(lambda row: calculate_average_stoichiometry(row, cond), axis=1)
        df['success' + cond] = df['coverage_' + cond] * df['stoichiometry_' + cond]
        df['success' + cond] = df['success' + cond].apply(lambda x: np.round(x).tolist())
        df['failure' + cond] = df['coverage_' + cond] - df['success' + cond]
    df.reset_index(inplace=True, drop=True)
    # Select columns starting with 'avg'
    avg_columns = df.filter(regex='^avg_stoichi', axis=1)

    # Calculate the absolute difference and store it in a new column
    df['delta'] = avg_columns.diff(axis=1).abs().iloc[:, -1]
    return df


def apply_GLM(df):
    conditions = ["_".join(x.split("_")[1:]) for x in df.columns if x.startswith('coverage_')]
    for index, row in df.iterrows():
        # print("processing index",index)
        coverage_columns = [col for col in row.index if col.startswith('coverage_')]
        conditions_reps = [c  for c in conditions for item in row[f'stoichiometry_{c}']]
        if(row[coverage_columns].apply(one_non_zero).sum()==len(coverage_columns)):
            # print(row)
            coverage = [item for c in conditions for item in row[f'coverage_{c}']]
            stoichiometry = [item  for c in conditions for item in row[f'stoichiometry_{c}']]
            #conditions_reps = [c  for c in conditions for item in row[f'stoichiometry_{c}']]
            transformed_row = pd.DataFrame({'condition': conditions_reps,'coverage': coverage,'stoichiometry': stoichiometry})
            pvalue = test_group_effect_bernoulli(transformed_row)
        else:
            fstat, pvalue = test_group_effect_binom(row, conditions_reps)

        df.at[index, 'pval'] = pvalue
        if row['delta'] < args.Delta_threshold:
           df.at[index, 'adj_pval'] = 1
        else:
           df.at[index, 'adj_pval'] = 0
    return df




if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--data_file',
                        help='tab-separated text files with first column M2_file_path, second column RepName, third column Condition')
    parser.add_argument('-o','--output_file', help='Output table.')
    parser.add_argument('-p', '--M2_threshold', type=float, default=0,
                        help='Compare only sites with at least one sample modified with site-probability more than p. Default =0')
    parser.add_argument('-s', '--Stoich_threshold', type=float, default=0,
            help='Compare only sites with at least one sample modified with stoichiometry more than s. Default =0')
    parser.add_argument('-delta', '--Delta_threshold', type=float, default=0,
                        help='Adjust only for sites with at least two conditions with stoichiometry difference greater than delta. Default =0.0')
    parser.add_argument('--ncpus', '-n', help='Number of CPU threads', type=int, default=1)
    args = parser.parse_args()
    warnings.filterwarnings("ignore")

    df_input = combine_samples(args)
    ncpus = args.ncpus
    df_list = np.array_split(df_input, ncpus)
    p = Pool(ncpus)
    time.sleep(1)
    pval_df_list = p.map(apply_GLM, df_list)
    df = pd.concat(pval_df_list)
    pvalues = df.loc[df['adj_pval'] == 0, 'pval']
    rej, pval_corr = smm.multipletests(pvalues, alpha=0.05, method='fdr_bh')[:2]
    df.loc[df['adj_pval'] == 0, 'adj_pval'] = pval_corr
    columns_to_drop = df.filter(regex='^success|^failure|^delta', axis=1).columns.tolist()
    df.drop(columns_to_drop, axis=1, inplace=True)
    df.to_csv(args.output_file, sep="\t", index=False)
    time.sleep(1)
    p.close()


