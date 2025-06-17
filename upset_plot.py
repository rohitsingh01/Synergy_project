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

synergy_responders = set(['A375', 'A498', 'A549', 'COLO741', 'Hs294T', 'MCF10A', 'MDAMB468',
 'MDST8', 'NCIH2170', 'NUGC3', 'OVISE', 'Panc0403', 'SBC5', 'THP1', 'YD10B']
)

synergy_nonresponders = set(['A2058', 'ACHN', 'AsPC1', 'C3A', 'CAL27', 'COLO201', 'COV434', 'ChaGoK1',
 'G401', 'HARA', 'HCC15', 'HCT116', 'HUTU80', 'HeLa', 'HuH1', 'HuH7',
 'KMS26', 'LMSU', 'LS411N', 'Li7', 'MKN45', 'PANC1', 'SCC25', 'SKHEP1',
 'SKMEL24', 'SKUT1', 'SUDHL10', 'SUDHL2', 'SUDHL4', 'SUIT2',
 'SUM159PT', 'SW1573', 'SW620', 'TE11', 'TE5', 'WILL1']
)

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