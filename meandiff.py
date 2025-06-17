import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load your data
df = pd.read_csv("all_gene_differences.csv")  # or load manually if needed

# Drop rows with missing PValues
df = df.dropna(subset=["PValue"])

# Compute -log10(p-value)
df["neg_log10_p"] = -np.log10(df["PValue"])

# Plot
plt.figure(figsize=(10, 7))
plt.scatter(df["MeanDiff"], df["neg_log10_p"], color='gray', alpha=0.7)

# Highlight significant points
sig = (df["PValue"] < 0.05) & (abs(df["MeanDiff"]) > 0.1)
plt.scatter(df[sig]["MeanDiff"], df[sig]["neg_log10_p"], color='red', label='Significant')

plt.axhline(-np.log10(0.05), color='blue', linestyle='--', label='p=0.05')
plt.axvline(0.1, color='green', linestyle='--')
plt.axvline(-0.1, color='green', linestyle='--')

plt.xlabel("Mean Difference")
plt.ylabel("-log10(P-Value)")
plt.title("Volcano Plot of Gene Expression")
plt.legend()
plt.tight_layout()
plt.show()
