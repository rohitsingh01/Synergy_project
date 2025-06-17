import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import ttest_ind
from matplotlib.patches import Patch
from io import StringIO
from matplotlib.gridspec import GridSpec

def normalize(name):
    return name.upper().replace('-', '').replace(' ', '')


folder = 'GSEA_INTER_GENES_CSV'
csv_files = [f for f in os.listdir(folder) if f.endswith('.csv')]


cell_lines = [
    'A2058', 'A549', 'HEPG2', 'CAL-27', 'HARA', 'Hs 294T', 'HuH-7', 'SKLMS1', 'MDA-MB-468', 'MDST8',
    'NCI-H2170', 'OVISE', 'PANC-1', 'SK-HEP-1', 'SK-MEL-24', 'SK-UT-1', 'SUIT-2', 'TE-11', 'TE-5',
    'THP-1', 'A375', 'AsPC-1', 'COLO-201', 'COV434', 'ChaGo-K-1', 'G-401', 'HCC-15', 'HCT-116',
    'HuH-1', 'KMS-26', 'LMSU', 'LS-411N', 'Li-7', 'MKN45', 'NUGC-3', 'Panc0403', 'SCC-25',
    'SU-DHL-10', 'SU-DHL-4', 'SUM159PT', 'SW1573', 'SW620', 'WILL-1'
]
normalized_cell_lines = [normalize(x) for x in cell_lines]
expr_df = pd.DataFrame({'CellLine': normalized_cell_lines})

# Merge all gene expression CSVs into one DataFrame
for file in csv_files:
    df = pd.read_csv(os.path.join(folder, file))
    gene_name = file.split()[0]
    df['CellLine'] = df['Cell Line Name'].apply(normalize)
    df_subset = df[['CellLine', 'Expression Public 24Q4']].rename(columns={'Expression Public 24Q4': gene_name})
    expr_df = pd.merge(expr_df, df_subset, on='CellLine', how='left')

expr_df.set_index('CellLine', inplace=True)
expr_log = np.log1p(expr_df)

# Define responders
responders = set(normalize(x) for x in [
    'A2058', 'A549', 'HEPG2', 'CAL-27', 'HARA', 'Hs 294T', 'HuH-7', 'SKLMS1', 'MDA-MB-468',
    'MDST8', 'NCI-H2170', 'OVISE', 'PANC-1', 'SK-HEP-1', 'SK-MEL-24', 'SK-UT-1', 'SUIT-2',
    'TE-11', 'TE-5', 'THP-1'
])

# Build metadata
meta = pd.DataFrame(index=expr_log.index)
meta['Response'] = ['Responder' if idx in responders else 'Non-responder' for idx in meta.index]

# Add cancer types
cancer_data = """
A549,Lung
HCT-116,Colorectal
HEPG2,Liver
CAL-27,Head and Neck
HARA,Lung
Hs 294T,Skin
HuH-7,Liver
SKLMS1,Soft tissue
MDA-MB-468,Breast
MDST8,Colorectal
NCI-H2170,Lung
OVISE,Ovary
PANC-1,Pancreas
SK-HEP-1,Liver
SK-MEL-24,Skin
SK-UT-1,Soft tissue
SUIT-2,Pancreas
TE-11,Esophagus
TE-5,Esophagus
THP-1,Leukemia
A375,Skin
AsPC-1,Pancreas
COLO-201,Colorectal
COV434,Ovary
ChaGo-K-1,Lung
G-401,Kidney
HCC-15,Lung
HuH-1,Liver
KMS-26,Plasma cell myeloma
LMSU,Gastric
LS-411N,Colorectal
Li-7,Liver
MKN45,Gastric
NUGC-3,Gastric
Panc0403,Pancreas
SCC-25,Head and Neck
SU-DHL-10,Lymphoma
SU-DHL-4,Lymphoma
SUM159PT,Breast
SW1573,Lung
SW620,Colorectal
WILL-1,Lymphoma
A2058,Skin
"""
cancer_df = pd.read_csv(StringIO(cancer_data), header=None, names=["CellLine", "CancerType"])
cancer_df["CellLine"] = cancer_df["CellLine"].apply(normalize)
cancer_df.set_index("CellLine", inplace=True)
meta["CancerType"] = cancer_df["CancerType"]


scaler = MinMaxScaler()
expr_scaled = pd.DataFrame(
    scaler.fit_transform(expr_log),
    index=expr_log.index,
    columns=expr_log.columns
)

# Compute differences between responders and non-responders
responder_idx = meta[meta['Response'] == 'Responder'].index
nonresponder_idx = meta[meta['Response'] == 'Non-responder'].index
t_stat, p_vals = ttest_ind(expr_scaled.loc[responder_idx], expr_scaled.loc[nonresponder_idx], axis=0, equal_var=False)

diff_df = pd.DataFrame({
    'Gene': expr_scaled.columns,
    'MeanDiff': expr_scaled.loc[responder_idx].mean() - expr_scaled.loc[nonresponder_idx].mean(),
    'PValue': p_vals
}).set_index('Gene')

top_genes = diff_df.sort_values(by='PValue').head(100).index
top_expr = expr_scaled[top_genes]


meta['GroupOrder'] = meta['Response'].map({'Non-responder': 0, 'Responder': 1})
meta['MeanExpr'] = top_expr.mean(axis=1)
sorted_idx = meta.sort_values(by=['GroupOrder', 'MeanExpr'], ascending=[True, False]).index
top_expr_sorted = top_expr.loc[sorted_idx]


response_palette = {'Responder': '#377EB8', 'Non-responder': '#E41A1C'}
cancer_palette = dict(zip(meta["CancerType"].unique(), sns.color_palette("hls", len(meta["CancerType"].unique()))))

col_colors = pd.DataFrame({
    'Response': meta["Response"].map(response_palette),
    'Cancer Type': meta["CancerType"].map(cancer_palette)
}, index=meta.index)
col_colors_sorted = col_colors.loc[sorted_idx]

sns.set(style="white")
g = sns.clustermap(
    top_expr_sorted.T,
    col_colors=col_colors_sorted,
    col_cluster=False,
    row_cluster=False,
    cmap="Reds",
    figsize=(50, 40)
)


response_legend = [Patch(color=color, label=label) for label, color in response_palette.items()]
cancer_legend = [Patch(color=color, label=label) for label, color in cancer_palette.items()]
g.ax_heatmap.legend(
    handles=response_legend + cancer_legend,
    title='Cell Line Info',
    bbox_to_anchor=(1.1, 1),
    loc='upper left'
)
plt.suptitle("Non-Responders Show Higher Gene Expression", fontsize=18, y=1.05)
plt.xlabel("Cell Lines")
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)
plt.ylabel("Expression differential")
plt.tight_layout()
plt.subplots_adjust(left=0.1, right=0.7, top=0.9, bottom=0.2)
plt.show()

pd.set_option('display.max_rows', None)
ranked_diff = diff_df.sort_values(by='MeanDiff', ascending=False)
print(ranked_diff[['MeanDiff', 'PValue']])
print(len(ranked_diff))
ranked_diff.to_csv('all_gene_differences.csv')
