import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
#DDX19A Data
group_a = [
    6.105624549, 3.820278465, 4.996271849, 6.202028179, 5.728481994,
    6.192334311, 5.859422186, 6.157016965, 5.137821465, 5.305939832,
    5.674404309, 4.87452353, 5.045326777, 5.679542225, 4.704071324,
    4.471558863, 4.512104045, 4.626694846, 4.325737653, 4.589286023,
    5.481975479, 6.795158678, 5.353232939, 0.382706767
]
group_b = [
    5.621172753, 5.14689956, 5.854245054, 4.705425039, 5.709842019,
    5.550900665, 5.683977107, 4.91981677, 4.900625333, 5.172727518,
    5.503030646, 5.100977648, 4.861955364, 4.9800253, 5.550285049,
    5.64990282, 5.358607249, 4.497612366, 4.77452353, 5.116447936
]
df = pd.DataFrame({
    'Expression': group_a + group_b,
    'Group': ['Nonresponders'] * len(group_a) + ['Responders'] * len(group_b)
})
plt.figure(figsize=(7, 5))
sns.boxplot(data=df, x='Group', y='Expression', palette='Reds')
sns.stripplot(data=df, x='Group', y='Expression', color='black', size=4, jitter=True, alpha=0.6)

plt.title('Comparison of DDX19A Expression Levels')
plt.ylabel('Expression Level')
plt.grid(axis='y')
plt.tight_layout()
plt.show()

#MAVS data
Mavs_nonresponders =   [4.319555769, 4.621678952, 5.559431619, 3.706442228, 4.702290585,
    4.004002316, 3.289033824, 4.926802684, 3.57118746, 4.243230135,
    3.702884409, 4.041106311, 4.489566812, 3.900123353, 3.33572706,
    2.837686761, 2.987525271, 3.74038956, 3.702884409, 3.766756592,
    4.541616269, 4.559431619, 3.490942773, 0.407137995]
Mavs_responders =  [
    3.993674362, 3.67468662, 4.339137385, 4.099295204, 4.300855871,
    3.6622055, 4.552745937, 4.876271361, 4.247927513, 4.976821852,
    4.22881869, 4.176322773, 3.648465443, 3.371558863, 4.244125943,
    3.937344392, 3.337711092, 3.058316496, 3.214124805, 3.218781168
]

df2 = pd.DataFrame({
    'Expression': Mavs_nonresponders + Mavs_responders,
    'Group': ['Nonresponders'] * len(group_a) + ['Responders'] * len(group_b)
})
plt.figure(figsize=(7, 5))
sns.boxplot(data=df2, x='Group', y='Expression', palette='Reds')
sns.stripplot(data=df2, x='Group', y='Expression', color='black', size=4, jitter=True, alpha=0.6)

plt.title('Comparison of MAVS Expression Levels')
plt.ylabel('Expression Level')
plt.grid(axis='y')
plt.tight_layout()
plt.show()

#TBK1 boxplot


tbk1_nonresponders =    [5.285073638, 4.569234794, 4.970857864, 4.746162657, 4.305548911,
    5.574436043, 5.642567608, 5.763059924, 5.441630009, 4.841466986,
    4.94899843, 5.299672345, 5.416869805, 5.159614856, 5.052333566,
    4.973321063, 3.535628594, 4.329587923, 4.31800615, 4.674707046,
    5.975042803, 4.484740587, 5.165227623, 0.020677353]
tbk1_responders = [
    4.853496704, 5.463360886, 3.92314918, 5.04176865, 5.015693807,
    5.165510018, 4.453517579, 4.711494907, 4.276496666, 4.024142346,
    4.343407822, 3.94204526, 5.264911693, 4.796493929, 4.799605422,
    5.010779839, 4.960697039, 4.007195501, 4.22342255, 4.514753498
]

df2 = pd.DataFrame({
    'Expression': tbk1_nonresponders + tbk1_responders,
    'Group': ['Nonresponders'] * len(group_a) + ['Responders'] * len(group_b)
})
plt.figure(figsize=(7, 5))
sns.boxplot(data=df2, x='Group', y='Expression', palette='Reds')
sns.stripplot(data=df2, x='Group', y='Expression', color='black', size=4, jitter=True, alpha=0.6)

plt.title('Comparison of TBK1 Expression Levels')
plt.ylabel('Expression Level')
plt.grid(axis='y')
plt.tight_layout()
plt.show()


#Sting boxplots

sting_nonresponders = [
    4.839308071, 6.458255587, 3.826831217, 0.421928095, 1.873996325,
    4.472255648, 6.139357765, 2.917623258, 0.554175893, 2.13562391,
    4.760495132, 3.059770155, 3.562052319, 1.177242999, 3.421928095,
    3.32650853, 4.683759754, 0.301633861, 0.363034406, 0.971843649,
    1.6360529, 2.859155834, 3.629820947, 0.085715951
]
sting_responders = [
    4.347665656, 0.687060688, 1.521050737, 4.493775408, 4.368069877,
    4.585563498, 1.028569152, 6.114991613, 5.579542225, 2.198494154,
    2.769771739, 1.831877241, 3.308885057, 2.280956314, 5.360013042,
    1.769771739, 5.617944917, 3.909773104, 4.806324057, 7.988287039
]
df2 = pd.DataFrame({
    'Expression':  sting_nonresponders + sting_responders,
    'Group': ['Nonresponders'] * len(group_a) + ['Responders'] * len(group_b)
})
plt.figure(figsize=(7, 5))
sns.boxplot(data=df2, x='Group', y='Expression', palette='Reds')
sns.stripplot(data=df2, x='Group', y='Expression', color='black', size=4, jitter=True, alpha=0.6)

plt.title('Comparison of STING Expression Levels')
plt.ylabel('Expression Level')
plt.grid(axis='y')
plt.tight_layout()
plt.show()

#Rigi

rigi_responders = [
    0.787060688,
    4.162639828,
    4.270726276,
    2.052333566,
    0.713531653,
    2.360025656,
    2.618535139,
    2.59057013,
    1.782573297,
    0.45614381,
    1.873996325,
    1.197610797,
    1.878208576,
    1.450497247,
    3.261887682,
    3.346408087,
    3.401587647,
    1.211031312,
    0.741546029,
    1.50599236,
    1.276322773,
    1.24404637,
    2.766756592,
    0.202462257
]
rigi_nonresponders = [
    1.718087584,
    1.967168608,
    1.361768359,
    5.915281866,
    1.765534746,
    2.327687364,
    1.66448284,
    2.107687869,
    2.927896454,
    1.459431619,
    1.014355293,
    3.321928095,
    3.752748591,
    1.843983844,
    1.929790998,
    2.901108243,
    1.495695163,
    2.510961919,
    3.424922088,
    1.189033824
]

df2 = pd.DataFrame({
    'Expression':  rigi_nonresponders + rigi_responders,
    'Group': ['Nonresponders'] * len(group_a) + ['Responders'] * len(group_b)
})
plt.figure(figsize=(7, 5))
sns.boxplot(data=df2, x='Group', y='Expression', palette='Reds')
sns.stripplot(data=df2, x='Group', y='Expression', color='black', size=4, jitter=True, alpha=0.6)

plt.title('Comparison of RIGI Expression Levels')
plt.ylabel('Expression Level')
plt.grid(axis='y')
plt.tight_layout()
plt.show()
