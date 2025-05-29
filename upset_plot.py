import pandas as pd
from upsetplot import UpSet, from_memberships
import matplotlib.pyplot as plt
dadt_responders = set([
    'A2058', 'A498', 'A549', 'C3A', 'CAL-27', 'COLO-741', 'HARA', 'HUTU-80',
    'Hs 294T', 'HuH-7', 'KLM-1', 'MDA-MB-468', 'MDST8', 'NCI-H2170', 'OVISE',
    'PANC-1', 'SBC-5', 'SK-HEP-1', 'SK-MEL-24', 'SK-UT-1', 'SUIT-2', 'TE-11', 'TE-5',
    'THP-1', 'YD-10B'
])

dadt_nonresponders = set([
    'A375', 'ACHN', 'AsPC-1', 'COLO-201', 'COV434', 'ChaGo-K-1', 'G-401', 'HCC-15',
    'HCT-116', 'HeLa', 'HuH-1', 'KMS-26', 'LMSU', 'LS-411N', 'Li-7', 'MCF-10A',
    'MKN45', 'NUGC-3', 'Panc 04.03', 'SCC-25', 'SU-DHL-10', 'SU-DHL-2', 'SU-DHL-4',
    'SUM159PT', 'SW1573', 'SW620', 'WILL-1'
])

synergy_responders = set([
    'A375', 'A549', 'MDA-MB-468', 'MDST8', 'OVISE', 'Panc 04.03', 'SBC-5', 'THP-1', 'YD-10B'
])

synergy_nonresponders = set([
    'A2058', 'A498', 'ACHN', 'AsPC-1', 'C3A', 'CAL-27', 'COLO-201', 'COLO-741',
    'COV434', 'ChaGo-K-1', 'G-401', 'HARA', 'HCC-15', 'HCT-116', 'HUTU-80', 'HeLa',
    'Hs 294T', 'HuH-1', 'HuH-7', 'KLM-1', 'KMS-26', 'LMSU', 'LS-411N', 'Li-7',
    'MCF-10A', 'MKN45', 'NCI-H2170', 'NUGC-3', 'PANC-1', 'SCC-25', 'SK-HEP-1',
    'SK-MEL-24', 'SK-UT-1', 'SU-DHL-10', 'SU-DHL-2', 'SU-DHL-4', 'SUIT-2',
    'SUM159PT', 'SW1573', 'SW620', 'TE-11', 'TE-5', 'WILL-1'
])

all_cell_lines = sorted(dadt_responders | dadt_nonresponders | synergy_responders | synergy_nonresponders)
data = []
for cell in all_cell_lines:
    membership = (
        ('DADT_Responder' if cell in dadt_responders else 'DADT_NonResponder'),
        ('Synergy_Responder' if cell in synergy_responders else 'Synergy_NonResponder')
    )
    data.append(membership)
upset_data = from_memberships(data)
plt.figure(figsize=(12, 6))
UpSet(upset_data, subset_size='count', show_counts=True).plot()
plt.title('UpSet Plot: DADT and Synergy Response Intersections')
plt.show()