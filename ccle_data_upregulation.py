import pandas as pd
import numpy as np

# Load and preprocess data
df = pd.read_csv('CCLE_data.csv')
df.drop(columns=['Name'], inplace=True)
df.set_index('Description', inplace=True)

# List of responder cell lines
responders = [
    'A2058', 'A498', 'A549', 'C3A', 'COLO-741', 'HARA', 'HUTU-80',
    'Hs 294T', 'KLM-1', 'MDA-MB-468', 'MDST8', 'NCI-H2170', 'OVISE',
    'PANC-1', 'SBC-5', 'SK-HEP-1', 'SK-UT-1', 'SUIT-2', 'TE-11',
    'TE-5', 'THP-1', 'YD-10B', 'CAL-27', 'HuH-7', 'SK-MEL-24'
]


responder_cols = [col for col in df.columns if any(r in col for r in responders)]

df_log = np.log2(df + 1)
responder_means = df_log[responder_cols].mean(axis=1)
threshold = responder_means.median()

upregulated_genes = responder_means[responder_means > threshold].sort_values(ascending=False)
downregulated_genes = responder_means[responder_means< threshold].sort_values(ascending=False)
print("Top upregulated genes in responders:")
print(upregulated_genes.head(10))
print(f"\nTotal upregulated genes in responders: {len(upregulated_genes)}")

upregulated_genes.to_csv('upregulated_in_responders.csv')
downregulated_genes.to_csv('downregulated.csv')