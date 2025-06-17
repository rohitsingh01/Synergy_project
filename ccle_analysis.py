import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load your CSV, make sure 'Name' is index or at least unique identifier
df = pd.read_csv("CCLE_data.csv")

# Set 'Description' or 'Name' as index (gene identifier)
df.set_index("Description", inplace=True)

# Select cell lines / tissue columns of interest (example list based on your columns)
selected_cell_lines = [
    "22RV1_PROSTATE", "2313287_STOMACH", "253JBV_URINARY_TRACT", "253J_URINARY_TRACT",
    "42MGBA_CENTRAL_NERVOUS_SYSTEM", "5637_URINARY_TRACT", "59M_OVARY",
    "639V_URINARY_TRACT", "647V_URINARY_TRACT", "697_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
    "769P_KIDNEY", "786O_KIDNEY", "8305C_THYROID", "8505C_THYROID"
]

# Filter the dataframe to include only those columns
filtered_df = df[selected_cell_lines]

# Optional: Check if all columns exist
missing_cols = [col for col in selected_cell_lines if col not in df.columns]
if missing_cols:
    print("Warning: Missing columns in data:", missing_cols)

# Plot heatmap of expression values
plt.figure(figsize=(9, max(6, len(filtered_df) * 0.2)))  # size tuned for many genes

sns.heatmap(filtered_df, cmap="coolwarm", annot=False,
            xticklabels=True, yticklabels=False)  # yticklabels off if too many genes

plt.title("Gene Expression Heatmap Across Selected Cell Lines")
plt.xlabel("Cell Lines / Tissues")
plt.ylabel("Genes")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
