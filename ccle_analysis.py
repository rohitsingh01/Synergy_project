import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap.umap_ as umap
df = pd.read_csv("CCLE_data.csv", index_col=0)
def normalize(name):
    return name.upper().replace("-", "").replace(" ", "")

nonresponders = [
    "A375", "ASPC1", "COLO201", "COV434", "CHAGOK1", "G401", "HCC15",
    "HCT116", "HUH1", "KMS26", "LMSU", "LS411N", "LI7", "MKN45", "NUGC3",
    "PANC0403", "SCC25", "SUDHL10", "SUDHL4", "SUM159PT", "SW1573", "SW620", "WILL1"
]
responders = [
    'A2058', 'A498', 'A549', 'C3A', 'COLO-741', 'HARA', 'HUTU-80',
    'Hs 294T', 'KLM-1', 'MDA-MB-468', 'MDST8', 'NCI-H2170', 'OVISE',
    'PANC-1', 'SBC-5', 'SK-HEP-1', 'SK-UT-1', 'SUIT-2', 'TE-11',
    'TE-5', 'THP-1', 'YD-10B','CAL-27', 'HuH-7', 'SK-MEL-24'
]

nonresponders_norm = set(normalize(x) for x in nonresponders)
responders_norm = set(normalize(x) for x in responders)

col_map = {col: normalize(col.split("_")[0]) for col in df.columns}
selected_cols = [col for col in df.columns if col_map.get(col) in nonresponders_norm.union(responders_norm)]
df_selected = df[selected_cols]
gene_var = df_selected.var(axis=1)
top_genes = gene_var.sort_values(ascending=False).head(100).index
df_top = df_selected.loc[top_genes]
df_log = np.log1p(df_top)
X = df_log.T.values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

labels = [("Nonresponder" if col_map[col] in nonresponders_norm else "Responder") for col in df_log.columns]


pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# plt.figure(figsize=(8,6))
# sns.scatterplot(x=X_pca[:,0], y=X_pca[:,1], hue=labels, palette={"Responder":"blue", "Nonresponder":"red"})
# plt.title("PCA (Standardized) of Top 500 Variable Genes")
# plt.xlabel("PC1")
# plt.ylabel("PC2")
# plt.legend()
# plt.show()

# UMAP
reducer = umap.UMAP(
    n_neighbors=5,
    min_dist=0.1,
    metric='correlation',
    random_state=42
)
X_umap = reducer.fit_transform(X_scaled)


# plt.figure(figsize=(8,6))
# sns.scatterplot(x=X_umap[:,0], y=X_umap[:,1], hue=labels, palette={"Responder":"blue", "Nonresponder":"red"})
# plt.title("UMAP (Standardized) of Top 500 Variable Genes")
# plt.xlabel("UMAP1")
# plt.ylabel("UMAP2")
# plt.legend()
# plt.show()

# from sklearn.metrics import silhouette_score

numeric_labels = [0 if l == "Nonresponder" else 1 for l in labels]

# score_pca = silhouette_score(X_pca, numeric_labels)
# score_umap = silhouette_score(X_umap, numeric_labels)
#
# print(f"Silhouette Score (PCA): {score_pca:.3f}")
# print(f"Silhouette Score (UMAP): {score_umap:.3f}")
# from sklearn.manifold import TSNE
#
# tsne = TSNE(n_components=2, perplexity=10, metric='cosine', random_state=42)
# X_tsne = tsne.fit_transform(X_scaled)
#
from sklearn.metrics import silhouette_score
# score_tsne = silhouette_score(X_tsne, numeric_labels)
# print(f"Silhouette Score (t-SNE): {score_tsne:.3f}")
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
lda = LDA(n_components=1)
X_lda = lda.fit_transform(X_scaled, numeric_labels)

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Create DataFrame with LDA1 and label
lda_df = pd.DataFrame({
    "LDA1": X_lda[:, 0],
    "Label": labels
})

lda_df["y"] = 0

lda_df = pd.DataFrame({
    "LDA1": X_lda[:, 0],
    "Label": labels
})
lda_df["y"] = 0  # dummy y

# Calculate group means and 95% CI
group_means = lda_df.groupby("Label")["LDA1"].mean()
# Plot
plt.figure(figsize=(10, 2.5))
sns.stripplot(data=lda_df, x="LDA1", y="y", hue="Label", jitter=0.2, dodge=True,
              palette={"Responder": "blue", "Nonresponder": "red"}, alpha=0.7)

# Draw mean lines and shaded CIs
# for label, color in zip(["Responder", "Nonresponder"], ["blue", "red"]):
#     mean = group_means[label]
#     plt.axvline(mean, color=color, linestyle="--", linewidth=1.5, label=f"{label} Mean")
#     plt.fill_betweenx([-0.1, 0.1], color=color, alpha=0.2)

# Style
plt.title("LDA Scatterplot with Group Means and 95% Confidence Intervals")
plt.xlabel("LDA Axis 1 (Discriminant Axis)")
plt.xticks(fontsize=7, rotation = 90)

plt.legend(loc="upper right")
plt.grid(True, axis="x", linestyle="--", alpha=0.3)
plt.tight_layout()
plt.show()

lda_values = X_lda[:, 0]
print(lda_values)
responder_vals = lda_values[np.array(labels) == "Responder"]
nonresponder_vals = lda_values[np.array(labels) == "Nonresponder"]
from scipy.stats import ttest_ind

t_stat, p_val = ttest_ind(responder_vals, nonresponder_vals)
print(f"T-test: t = {t_stat:.3f}, p = {p_val:.4e}")
def cohens_d(a, b):
    mean_diff = np.mean(a) - np.mean(b)
    pooled_sd = np.sqrt((np.std(a, ddof=1)**2 + np.std(b, ddof=1)**2) / 2)
    return mean_diff / pooled_sd

d = cohens_d(responder_vals, nonresponder_vals)
print(f"Cohen's d = {d:.3f}")

#repeat with some random cell line
# look at upregulated responders
from scipy.stats import ttest_ind

responder_cols = [col for col in df_log.columns if col_map[col] in responders_norm]
nonresponder_cols = [col for col in df_log.columns if col_map[col] in nonresponders_norm]
responder_data = df_log[responder_cols]
nonresponder_data = df_log[nonresponder_cols]

results = []

for gene in df_log.index:
    r_vals = responder_data.loc[gene]
    nr_vals = nonresponder_data.loc[gene]
    t_stat, p_val = ttest_ind(r_vals, nr_vals, equal_var=False)
    fc = np.mean(r_vals) - np.mean(nr_vals)

    results.append((gene, fc, p_val))
results_df = pd.DataFrame(results, columns=["Gene", "log2FC", "p_value"])
from statsmodels.stats.multitest import multipletests

results_df["adj_p"] = multipletests(results_df["p_value"], method="fdr_bh")[1]
upregulated = results_df[(results_df["log2FC"] > 1) ]

print("Upregulated genes in responders:")
print(upregulated.sort_values("log2FC", ascending=False))

