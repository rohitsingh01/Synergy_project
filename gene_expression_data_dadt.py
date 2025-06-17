import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


folder_path = 'GENE EXPRESSION CSVs'
csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]


cell_lines = [
    'A2058', 'A549', 'HEPG2', 'CAL-27', 'HARA',
    'Hs 294T', 'HuH-7', 'SKLMS1', 'MDA-MB-468', 'MDST8', 'NCI-H2170', 'OVISE',
    'PANC-1', 'SK-HEP-1', 'SK-MEL-24', 'SK-UT-1', 'SUIT-2', 'TE-11', 'TE-5',
    'THP-1', 'A375', 'AsPC-1', 'COLO-201', 'COV434', 'ChaGo-K-1', 'G-401', 'HCC-15',
    'HCT-116', 'HuH-1', 'KMS-26', 'LMSU', 'LS-411N', 'Li-7',
    'MKN45', 'NUGC-3', 'Panc0403', 'SCC-25', 'SU-DHL-10', 'SU-DHL-4',
    'SUM159PT', 'SW1573', 'SW620', 'WILL-1'
]
def normalize(name):
    return name.upper().replace('-', '').replace(' ', '')

cell_lines_normalized = [normalize(cl) for cl in cell_lines]
expression_matrix = pd.DataFrame({'Cell Line Name': cell_lines_normalized})

for file in csv_files:
    df = pd.read_csv(os.path.join(folder_path, file))
    if 'Cell Line Name' not in df.columns or 'Expression Public 24Q4' not in df.columns:
        continue
    df['Cell Line Name'] = df['Cell Line Name'].apply(normalize)
    gene_name = file.split()[0]
    df_gene = df[['Cell Line Name', 'Expression Public 24Q4']].rename(columns={'Expression Public 24Q4': gene_name})
    expression_matrix = pd.merge(expression_matrix, df_gene, on='Cell Line Name', how='left')

expression_matrix.set_index('Cell Line Name', inplace=True)

# Define responders and non-responders
dadt_responders = set(map(normalize, [
    'A2058', 'A549', 'HEPG2', 'CAL-27', 'HARA',
    'Hs 294T', 'HuH-7', 'SKLMS1', 'MDA-MB-468', 'MDST8', 'NCI-H2170', 'OVISE',
    'PANC-1', 'SK-HEP-1', 'SK-MEL-24', 'SK-UT-1', 'SUIT-2', 'TE-11', 'TE-5',
    'THP-1',
]))
dadt_nonresponders = set(map(normalize, [
    'A375', 'AsPC-1', 'COLO-201', 'COV434', 'ChaGo-K-1', 'G-401', 'HCC-15',
    'HCT-116', 'HuH-1', 'KMS-26', 'LMSU', 'LS-411N', 'Li-7',
    'MKN45', 'NUGC-3', 'Panc0403', 'SCC-25', 'SU-DHL-10', 'SU-DHL-4',
    'SUM159PT', 'SW1573', 'SW620', 'WILL-1'
]))
print(len(dadt_nonresponders) + len(dadt_responders))

responders_matrix = expression_matrix.loc[expression_matrix.index.intersection(dadt_responders)]
nonresponders_matrix = expression_matrix.loc[expression_matrix.index.intersection(dadt_nonresponders)]
nonresponders_matrix+=0.1
combined = pd.concat([responders_matrix, nonresponders_matrix])
vmin = combined.min().min()
vmax = combined.max().max()
plt.figure(figsize=(10, 16))
sns.heatmap(responders_matrix, cmap='Reds', linewidths=0.3,
            cbar_kws={'label': 'Expression'}, yticklabels=True,
            vmin=vmin, vmax=vmax, annot=False, fmt=".2f", annot_kws={"size": 7})
plt.title('DADT Responders')
plt.ylabel('Cell Line Name')
plt.xlabel('Gene')
plt.tick_params(axis='y', labelsize=10)
plt.tight_layout()
plt.show()

# Non-responders heatmap
plt.figure(figsize=(10, 16))
sns.heatmap(nonresponders_matrix, cmap='Reds', linewidths=0.3,
            cbar_kws={'label': 'Expression'}, yticklabels=True,
            vmin=vmin, vmax=vmax, annot=False, fmt=".2f", annot_kws={"size": 7})
plt.title('DADT Non-Responders')
plt.ylabel('Cell Line Name')
plt.xlabel('Gene')
plt.tick_params(axis='y', labelsize=10)
plt.tight_layout()
plt.show()

common_genes = responders_matrix.columns.intersection(nonresponders_matrix.columns)
responders_data = responders_matrix[common_genes]
nonresponders_data = nonresponders_matrix[common_genes]

# Updated t-test loop using TWO-tailed test, and dropping NaNs and 0s
p_values = {}
for gene in common_genes:
    r = responders_data[gene].dropna()
    n = nonresponders_data[gene].dropna()

    # Drop zero values to reduce noise
    r = r[r != 0].astype(float)
    n = n[n != 0].astype(float)

    if len(r) > 2 and len(n) > 2:
        stat, p_two_tailed = ttest_ind(r, n, equal_var=False)
        p_values[gene] = p_two_tailed  # Use full two-tailed p-value
    else:
        p_values[gene] = pd.NA

pval_df = pd.DataFrame.from_dict(p_values, orient='index', columns=['p-value'])
pval_df.dropna(inplace=True)

# Multiple testing correction
_, fdr_corrected, _, _ = multipletests(pval_df['p-value'].values, method='fdr_bh')
pval_df['FDR'] = fdr_corrected

# Sort and filter
pval_df.sort_values('p-value', inplace=True)
significant_genes = pval_df[pval_df['FDR'] < 0.05]

# Print results
print(pval_df)
print(significant_genes)
