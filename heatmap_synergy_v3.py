import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import ttest_ind
from matplotlib.patches import Patch
from io import StringIO

# Normalize function for consistent cell line naming
def normalize(name):
    return name.upper().replace('-', '').replace(' ', '')

folder = 'GENE EXPRESSION CSVs'
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

# Merge gene expression CSVs
for file in csv_files:
    df = pd.read_csv(os.path.join(folder, file))
    gene_name = file.split()[0]
    df['CellLine'] = df['Cell Line Name'].apply(normalize)
    df_subset = df[['CellLine', 'Expression Public 24Q4']].rename(columns={'Expression Public 24Q4': gene_name})
    expr_df = pd.merge(expr_df, df_subset, on='CellLine', how='left')

expr_df.set_index('CellLine', inplace=True)

# Define responders
responders = set(normalize(x) for x in [
    'A2058', 'A498', 'A549', 'C3A', 'COLO-741', 'HARA', 'HUTU-80',
    'Hs 294T', 'KLM-1', 'MDA-MB-468', 'MDST8', 'NCI-H2170', 'OVISE',
    'PANC-1', 'SBC-5', 'SK-HEP-1', 'SK-UT-1', 'SUIT-2', 'TE-11',
    'TE-5', 'THP-1', 'YD-10B', 'A375', 'HCC-15', 'MCF-10A', 'NUGC-3', 'Panc 04.03'
])


expr_log = np.log1p(expr_df)


meta = pd.DataFrame(index=expr_log.index)
meta['Response'] = ['Responder' if idx in responders else 'Non-responder' for idx in meta.index]


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

# MinMax scale
scaler = MinMaxScaler()
expr_scaled = pd.DataFrame(
    scaler.fit_transform(expr_log),
    index=expr_log.index,
    columns=expr_log.columns
)


responder_idx = meta[meta['Response'] == 'Responder'].index
nonresponder_idx = meta[meta['Response'] == 'Non-responder'].index
t_stat, p_vals = ttest_ind(expr_scaled.loc[responder_idx], expr_scaled.loc[nonresponder_idx], axis=0, equal_var=False)

diff_df = pd.DataFrame({
    'Gene': expr_scaled.columns,
    'MeanDiff': expr_scaled.loc[responder_idx].mean() - expr_scaled.loc[nonresponder_idx].mean(),
    'PValue': p_vals
}).set_index('Gene')


top_genes = diff_df.sort_values(by='PValue').index
top_expr = expr_scaled[top_genes]


custom_order = [
    'MDST8', 'THP-1', 'MDA-MB-468', 'A549', 'OVISE', 'A375', 'NUGC-3', 'Hs 294T', 'SK-UT-1',
    'NCI-H2170', 'TE-11', 'SUIT-2', 'HCC-15', 'A2058', 'TE-5', 'PANC-1', 'SK-HEP-1', 'HARA',
    'COV434', 'HuH-7', 'CAL-27', 'HuH-1', 'COLO-201', 'SW620', 'HCT-116', 'MKN45', 'SCC-25',
    'WILL-1', 'SK-MEL-24', 'AsPC-1', 'G-401', 'LMSU', 'SW1573', 'LS-411N', 'SU-DHL-4', 'SU-DHL-10',
    'ChaGo-K-1', 'KMS-26', 'Li-7', 'SUM159PT', 'HEPG2', 'SKLMS1', 'Panc0403'
]
custom_order_normalized = [normalize(x) for x in custom_order]
sorted_idx = [idx for idx in custom_order_normalized if idx in expr_scaled.index]


top_expr_sorted = top_expr.loc[sorted_idx].fillna(0)


col_colors = pd.DataFrame({
    'Response': meta["Response"].map({'Responder': '#377EB8', 'Non-responder': '#E41A1C'}),
    'Cancer Type': meta["CancerType"].map(
        dict(zip(meta["CancerType"].unique(), sns.color_palette("hls", len(meta["CancerType"].unique()))))
    )
}, index=meta.index)
col_colors_sorted = col_colors.loc[sorted_idx]


def label_cellline(idx):
    return f"{meta.loc[idx, 'Response'][:3].upper()}: {idx}"
top_expr_sorted.index = [label_cellline(i) for i in top_expr_sorted.index]
col_colors_sorted.index = top_expr_sorted.index

sns.set(style="white")
g = sns.clustermap(
    top_expr_sorted.T,
    col_colors=col_colors_sorted,
    col_cluster=False,
    row_cluster=True,
    cmap="Reds",
    figsize=(50, 40)
)

n_nonresponders = sum(meta.loc[sorted_idx, 'Response'] == 'Non-responder')
plt.axvline(x=n_nonresponders - 0.5, color='black', linewidth=2, linestyle='--')

# Legends
response_legend = [Patch(color=color, label=label) for label, color in {'Responder': '#377EB8', 'Non-responder': '#E41A1C'}.items()]
cancer_legend = [
    Patch(color=color, label=label)
    for label, color in dict(zip(meta["CancerType"].unique(), sns.color_palette("hls", len(meta["CancerType"].unique())))).items()
]
g.ax_heatmap.legend(
    handles=response_legend + cancer_legend,
    title='Cell Line Info',
    bbox_to_anchor=(1.3, 1),
    loc='upper left'
)


plt.suptitle("Gene Expression Differences Between Non-Responders and Responders", fontsize=18, y=1.05)
plt.xlabel("Cell Lines (Labeled by Response)")
plt.ylabel("Top 100 Differential Genes")
plt.tight_layout()
plt.subplots_adjust(left=0.1, right=0.7, top=0.9, bottom=0.2)
plt.show()

# ranked_diff = diff_df.sort_values(by='MeanDiff', ascending=False)
# print(ranked_diff[['MeanDiff', 'PValue']])
# print(f"Total genes analyzed: {len(ranked_diff)}")
# ranked_diff.to_csv('original_gene_diff.csv')
#
# to_export = expr_scaled
#
# responders_df = to_export.loc[meta[meta['Response'] == 'Responder'].index]
# nonresponders_df = to_export.loc[meta[meta['Response'] == 'Non-responder'].index]
#
# # Export to CSV
# responders_df.to_csv('responder_expression.csv')
# nonresponders_df.to_csv('nonresponder_expression.csv')
#
# print("Exported responder_expression.csv and nonresponder_expression.csv successfully.")
# # Combine both responders and non-responders into one DataFrame in your custom order
combined_heatmap_data = top_expr.loc[sorted_idx].fillna(0)

# Add Response label for clarity
combined_heatmap_data.insert(0, 'Response', meta.loc[sorted_idx, 'Response'].values)

# Export to CSV
combined_heatmap_data.to_csv('heatmap_combined_expression.csv')

print("Exported heatmap_combined_expression.csv with both responders and non-responders in custom order.")
