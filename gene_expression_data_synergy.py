import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

folder_path = 'C:/Users/u1241114/PycharmProjects/SYNERGY_PROJECT/GENE EXPRESSION CSVs'
csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
cell_lines = [
    'A2058', 'A498', 'A549', 'HepG2', 'CAL27', 'COLO741', 'HARA', 'HUTU80',
    'Hs294T', 'HuH7', 'KLM1', 'MDAMB468', 'MDST8', 'NCIH2170', 'OVISE',
    'PANC1', 'SBC5', 'SKHEP1', 'SKMEL24', 'SKUT1', 'SUIT2', 'TE11', 'TE5',
    'THP1', 'YD10B', 'A375', 'ACHN', 'AsPC1', 'COLO201', 'COV434', 'ChaGoK1', 'G401', 'HCC15',
    'HCT116', 'HeLa', 'HuH1', 'KMS26', 'LMSU', 'LS411N', 'Li7', 'MCF10A',
    'MKN45', 'NUGC3', 'Panc0403', 'SCC25', 'SUDHL10', 'SUDHL2', 'SUDHL4',
    'SUM159PT', 'SW1573', 'SW620', 'WILL1'
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

synergy_responders = set(map(normalize, ['A375', 'A498', 'A549', 'COLO741', 'Hs294T', 'MCF10A', 'MDAMB468',
 'MDST8', 'NCIH2170', 'NUGC3', 'OVISE', 'Panc0403', 'SBC5', 'THP1', 'YD10B']))
synergy_nonresponders = set(map(normalize, ['A2058', 'ACHN', 'AsPC1', 'C3A', 'CAL27', 'COLO201', 'COV434', 'ChaGoK1',
 'G401', 'HARA', 'HCC15', 'HCT116', 'HUTU80', 'HeLa', 'HuH1', 'HuH7',
 'KMS26', 'LMSU', 'LS411N', 'Li7', 'MKN45', 'PANC1', 'SCC25', 'SKHEP1',
 'SKMEL24', 'SKUT1', 'SUDHL10', 'SUDHL2', 'SUDHL4', 'SUIT2',
 'SUM159PT', 'SW1573', 'SW620', 'TE11', 'TE5', 'WILL1']))
responders_matrix = expression_matrix.loc[expression_matrix.index.intersection(synergy_responders)]
nonresponders_matrix = expression_matrix.loc[expression_matrix.index.intersection(synergy_nonresponders)]

# Plot responders heatmap
plt.figure(figsize=(10, 16))
sns.heatmap(responders_matrix, cmap='magma', linewidths=0,
            cbar_kws={'label': 'Expression'}, yticklabels=True)
plt.title('Synergy Responders')
plt.ylabel('Cell Line Name')
plt.xlabel('Gene')
plt.tick_params(axis='y', labelsize=10)
plt.tight_layout()
plt.show()

# Plot non-responders heatmap
plt.figure(figsize=(10, 16))
sns.heatmap(nonresponders_matrix, cmap='magma', linewidths=0,
            cbar_kws={'label': 'Expression'}, yticklabels=True)
plt.title('Synergy Non-Responders')
plt.ylabel('Cell Line Name')
plt.xlabel('Gene')
plt.tick_params(axis='y', labelsize=10)
plt.tight_layout()
plt.show()