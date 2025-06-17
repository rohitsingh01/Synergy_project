#need to analyze basal expression levels of cyokines and TBK1 interactors in responders vs nonresponders
#need to add DDX19A to previous analysis
#look into TBK1_v1, v2


import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

folder_path = 'C:/Users/u1241114/PycharmProjects/SYNERGY_PROJECT/cytokine_csvs'

csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
print(len(csv_files))
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


dadt_responders = set(map(normalize, [
    'A2058', 'A498', 'A549', 'C3A', 'CAL-27', 'COLO-741', 'HARA', 'HUTU-80',
    'Hs 294T', 'HuH-7', 'KLM-1', 'MDA-MB-468', 'MDST8', 'NCI-H2170', 'OVISE',
    'PANC-1', 'SBC-5', 'SK-HEP-1', 'SK-MEL-24', 'SK-UT-1', 'SUIT-2', 'TE-11', 'TE-5',
    'THP-1', 'YD-10B'
]))
dadt_nonresponders = set(map(normalize, [
    'A375', 'ACHN', 'AsPC-1', 'COLO-201', 'COV434', 'ChaGo-K-1', 'G-401', 'HCC-15',
    'HCT-116', 'HeLa', 'HuH-1', 'KMS-26', 'LMSU', 'LS-411N', 'Li-7', 'MCF-10A',
    'MKN45', 'NUGC-3', 'Panc 04.03', 'SCC-25', 'SU-DHL-10', 'SU-DHL-2', 'SU-DHL-4',
    'SUM159PT', 'SW1573', 'SW620', 'WILL-1'
]))

responders_matrix = expression_matrix.loc[expression_matrix.index.intersection(dadt_responders)]
nonresponders_matrix = expression_matrix.loc[expression_matrix.index.intersection(dadt_nonresponders)]

plt.figure(figsize=(10, 16))
sns.heatmap(responders_matrix, cmap='magma', linewidths=0,
            cbar_kws={'label': 'Expression'}, yticklabels=True)
plt.title('DADT Responders')
plt.ylabel('Cell Line Name')
plt.xlabel('Gene')
plt.tick_params(axis='y', labelsize=10)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 16))
sns.heatmap(nonresponders_matrix, cmap='magma', linewidths=0,
            cbar_kws={'label': 'Expression'}, yticklabels=True)
plt.title('DADT Non-Responders')
plt.ylabel('Cell Line Name')
plt.xlabel('Gene')
plt.tick_params(axis='y', labelsize=10)
plt.tight_layout()
plt.show()
