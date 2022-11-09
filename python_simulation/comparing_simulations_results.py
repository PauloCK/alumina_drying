import numpy as np
import pandas as pd


def get_df_name(df):
    name =[x for x in globals() if globals()[x] is df][0]
    return name

def add_suffix_to_df_cols(df):
    suffix = get_df_name(df)[:-3]
    df.columns = [f'{col}_{suffix}' if col != 'Tempo_min' else f'{col}' for col in df.columns]
    return df

def create_comparison_cols(df, var, pos):
    df[f'diff_{var}_{pos}cm_p_s'] = 100*(abs(df[f'{var}_{pos}cm_python'] - df[f'{var}_{pos}cm_scilab'])) / np.maximum(abs(df[f'{var}_{pos}cm_python']), abs(df[f'{var}_{pos}cm_scilab']))
    df[f'diff_{var}_{pos}cm_p_m'] = 100*(abs(df[f'{var}_{pos}cm_python'] - df[f'{var}_{pos}cm_matlab'])) / np.maximum(abs(df[f'{var}_{pos}cm_python']), abs(df[f'{var}_{pos}cm_matlab']))
    df[f'diff_{var}_{pos}cm_s_m'] = 100*(abs(df[f'{var}_{pos}cm_scilab'] - df[f'{var}_{pos}cm_matlab'])) / np.maximum(abs(df[f'{var}_{pos}cm_scilab']), abs(df[f'{var}_{pos}cm_matlab']))

    return df

def create_statistics_df(df):
    suffix = df.columns[0][-4:]
    statistics_df = pd.DataFrame()

    for var in ['X', 'T']:
        for pos in [1,5,10]:
            var_pos = f'{var}_{pos}cm'
            statistics_df.at[var_pos, 'min_diff'] = min(df[f'diff_{var_pos}{suffix}'][2:])
            statistics_df.at[var_pos, 'max_diff'] = max(df[f'diff_{var_pos}{suffix}'])
            statistics_df.at[var_pos, 'mean_diff'] = np.mean(df[f'diff_{var_pos}{suffix}'])
            statistics_df.at[var_pos, 'median_diff'] = np.median(df[f'diff_{var_pos}{suffix}'])

    return statistics_df


scilab_df = pd.read_fwf('simulations_results/simulation_results_scilab.txt')
matlab_df = pd.read_fwf('simulations_results/simulation_results_matlab.txt')
python_df = pd.read_csv('simulations_results/simulation_results_python.csv')
dfs = [scilab_df, matlab_df, python_df]

for df in dfs:
    add_suffix_to_df_cols(df)

concat_df = pd.concat(dfs, axis=1)
for var in ['X', 'T']:
    for pos in [1,5,10]:
        create_comparison_cols(concat_df, var, pos)

diffs_s_m_df = concat_df[[col for col in concat_df.columns if 's_m' in col]]
diffs_p_m_df = concat_df[[col for col in concat_df.columns if 'p_m' in col]]
diffs_p_s_df = concat_df[[col for col in concat_df.columns if 'p_s' in col]]
diffs_dfs = [diffs_s_m_df, diffs_p_m_df, diffs_p_s_df]
for df in diffs_dfs:
    suffix = df.columns[0][-4:]
    statistics_df = create_statistics_df(df)
    statistics_df.to_csv(f'statistics_df{suffix}.csv')