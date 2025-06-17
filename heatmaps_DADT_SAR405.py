import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

df = pd.read_csv("HTRF_Normalised_Data_Full Screen.csv")

red_white = LinearSegmentedColormap.from_list('red_white', ['white', 'red'])
blue_white = LinearSegmentedColormap.from_list('blue_white', ['white', 'blue'])
purple_white = LinearSegmentedColormap.from_list('purple_white', ['white', 'purple'])

cell_lines_ec50_2 = df.groupby('cell_line')['ec50_2'].first().reset_index()
cell_lines_ec50_2 = cell_lines_ec50_2.sort_values(by='ec50_2', ascending=False)
ec50_2_top2 = cell_lines_ec50_2.iloc[:2]
ec50_2_rest = cell_lines_ec50_2.iloc[2:]
heatmap_data_2_rest = ec50_2_rest.set_index('cell_line').T
heatmap_data_2_top2 = ec50_2_top2.set_index('cell_line').T

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 3), gridspec_kw={'width_ratios': [len(heatmap_data_2_rest.columns), len(heatmap_data_2_top2.columns)]})
sns.heatmap(heatmap_data_2_rest, cmap=red_white, ax=ax1, annot=False, fmt=".2f", cbar_kws={'label': 'EC50'})
for idx, cell_line in enumerate(heatmap_data_2_rest.columns):
    val = heatmap_data_2_rest[cell_line].values[0]
    ax1.text(idx + 0.5, -0.3, f'{val:.2f}', ha='center', va='bottom', rotation=90, fontsize=8)
ax1.set_title('EC50 of DADT (Normal Range)')
ax1.set_ylabel('')
ax1.set_xlabel('Cell Line')
ax1.set_xticks(ax1.get_xticks())
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
ax1.set_ylim(1.2, -0.2)

sns.heatmap(heatmap_data_2_top2, cmap=red_white, ax=ax2, annot=False, fmt=".2f", cbar=False)
for idx, cell_line in enumerate(heatmap_data_2_top2.columns):
    val = heatmap_data_2_top2[cell_line].values[0]
    ax2.text(idx + 0.5, -0.3, f'{val:.2f}', ha='center', va='bottom', rotation=90, fontsize=8)
ax2.set_title('EC50 of DADT (Top 2 High)')
ax2.set_ylabel('')
ax2.set_xlabel('Cell Line')
ax2.set_xticks(ax2.get_xticks())
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90)
ax2.set_ylim(1.2, -0.2)

plt.tight_layout()
plt.show()

print(df.groupby('cell_line')['ec50_2'].first().reset_index().query('ec50_2 < 4')['cell_line'].values)
print(df.groupby('cell_line')['ec50_2'].first().reset_index().query('ec50_2 > 4')['cell_line'].values)

cell_lines_ec50_1 = df.groupby('cell_line')['ec50_1'].first().reset_index()
cell_lines_ec50_1 = cell_lines_ec50_1.sort_values(by='ec50_1', ascending=False)
ec50_1_top5 = cell_lines_ec50_1.iloc[:5]
ec50_1_rest = cell_lines_ec50_1.iloc[5:]
heatmap_data_1_rest = ec50_1_rest.set_index('cell_line').T
heatmap_data_1_top5 = ec50_1_top5.set_index('cell_line').T

fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(16, 3), gridspec_kw={'width_ratios': [len(heatmap_data_1_rest.columns), len(heatmap_data_1_top5.columns)]})
sns.heatmap(heatmap_data_1_rest, cmap=blue_white, ax=ax3, annot=False, fmt=".2f", cbar_kws={'label': 'EC50'})
for idx, cell_line in enumerate(heatmap_data_1_rest.columns):
    val = heatmap_data_1_rest[cell_line].values[0]
    ax3.text(idx + 0.5, -0.3, f'{val:.2f}', ha='center', va='bottom', rotation=90, fontsize=8)
ax3.set_title('EC50 of SAR405 (Normal Range)')
ax3.set_ylabel('')
ax3.set_xlabel('Cell Line')
ax3.set_xticks(ax3.get_xticks())
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=90)
ax3.set_ylim(1.2, -0.2)

sns.heatmap(heatmap_data_1_top5, cmap=blue_white, ax=ax4, annot=False, fmt=".2f", cbar=False)
for idx, cell_line in enumerate(heatmap_data_1_top5.columns):
    val = heatmap_data_1_top5[cell_line].values[0]
    ax4.text(idx + 0.5, -0.3, f'{val:.2f}', ha='center', va='bottom', rotation=90, fontsize=8)
ax4.set_title('EC50 of SAR405 (Top 5 High)')
ax4.set_ylabel('')
ax4.set_xlabel('Cell Line')
ax4.set_xticks(ax4.get_xticks())
ax4.set_xticklabels(ax4.get_xticklabels(), rotation=90)
ax4.set_ylim(1.2, -0.2)

plt.tight_layout()
plt.show()
cell_lines_synergy = df.groupby('cell_line')['Synergy_score'].first().reset_index()
cell_lines_synergy = cell_lines_synergy.sort_values(by='Synergy_score', ascending=False)

top_n = 5  # number of top synergy scores to separate
synergy_top = cell_lines_synergy.iloc[:top_n]
synergy_rest = cell_lines_synergy.iloc[top_n:]

heatmap_data_top = synergy_top.set_index('cell_line').T
heatmap_data_rest = synergy_rest.set_index('cell_line').T

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 3),
                               gridspec_kw={'width_ratios': [len(heatmap_data_rest.columns), len(heatmap_data_top.columns)]})

purple_white = LinearSegmentedColormap.from_list('purple_white', ['white', 'purple'])

# Plot rest (normal range)
sns.heatmap(heatmap_data_rest, cmap=purple_white, ax=ax1, cbar_kws={'label': 'Synergy Score'}, annot=False, fmt=".2f")
for idx, cell_line in enumerate(heatmap_data_rest.columns):
    val = heatmap_data_rest[cell_line].values[0]
    ax1.text(idx + 0.5, -0.3, f'{val:.2f}', ha='center', va='bottom', rotation=90, fontsize=8)
ax1.set_title('Synergy Score (Normal Range)')
ax1.set_ylabel('')
ax1.set_xlabel('Cell Line')
ax1.set_xticks(ax1.get_xticks())
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
ax1.set_ylim(1.2, -0.2)

# Plot top synergy scores
sns.heatmap(heatmap_data_top, cmap=purple_white, ax=ax2, cbar=False, annot=False, fmt=".2f")
for idx, cell_line in enumerate(heatmap_data_top.columns):
    val = heatmap_data_top[cell_line].values[0]
    ax2.text(idx + 0.5, -0.3, f'{val:.2f}', ha='center', va='bottom', rotation=90, fontsize=8)
ax2.set_title(f'Synergy Score (Top {top_n})')
ax2.set_ylabel('')
ax2.set_xlabel('Cell Line')
ax2.set_xticks(ax2.get_xticks())
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90)
ax2.set_ylim(1.2, -0.2)

plt.tight_layout()
plt.show()
print(df.groupby('cell_line')['Synergy_score'].first().reset_index().query('Synergy_score > 1')['cell_line'].values)
print(df.groupby('cell_line')['Synergy_score'].first().reset_index().query('Synergy_score < 1')['cell_line'].values)
