import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
#%%
input_path = "/Users/ireneu/berman_lab/Rett/Figure3_raw_results.csv"
data = pd.read_csv(input_path)
#%%
#3A
figure_data = data[data["figure"]=="3A"]

fig, ax = plt.subplots(figsize=(7,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group",
                 ci="sd",
                 errwidth=1)

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none'
                   )
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-2:], labels[-2:], loc=2, borderaxespad=0.)

box_pairs = [
    (('CD206', 'Con shRNA'), ('CD206', 'MeCP2 shRNA')),
    (('IL-13', 'Con shRNA'), ('IL-13', 'MeCP2 shRNA')),
    (('CD86', 'Con shRNA'), ('CD86', 'MeCP2 shRNA')),
    (('IL-1', 'Con shRNA'), ('IL-1', 'MeCP2 shRNA'))
]
s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch')

# test_results = add_stat_annotation(ax, data=figure_data, x='condition', y='result',
#                                    box_pairs=[('CD206', 'IL-13'),('CD86', 'IL-1')],
#                                    text_annot_custom=["M2 markers", "M1 markers"],
#                                    perform_stat_test=False, pvalues=[0, 0],
#                                    loc='outside', verbose=0)

plt.xlabel("")
plt.ylabel("Relative mRNA Expression")
plt.show()
# plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_3A.pdf', bbox_inches='tight')

#%%
figure_data = data[data["figure"]=="3B"]

fig, ax = plt.subplots(figsize=(3,5))
ax = sns.barplot(x='group', y='result', data=figure_data,
                 ci="sd",
                 errwidth=1)


# Add jitter with the swarmplot function
ax = sns.stripplot(x='group', y='result', data=figure_data, color="black",
                   jitter=True, marker='o', facecolors='none'
                   )


box_pairs = [
    ('Con shRNA', 'MeCP2 shRNA')
]

s = add_stat_annotation(ax, data=figure_data, x='group', y='result', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch')

plt.xlabel("")
plt.ylabel("% Phagocytosis")
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_3B.pdf', bbox_inches='tight')
#%%
#3C
figure_data = data[data["figure"]=="3C"]

fig, ax = plt.subplots(figsize=(4,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group",
                 ci="sd",
                 errwidth=1)

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none'
                   )
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-2:], labels[-2:], loc=1, borderaxespad=0.)

box_pairs = [
    (('GFAP', 'Con shRNA'), ('GFAP', 'MeCP2 shRNA')),
    (('MAP2', 'Con shRNA'), ('MAP2', 'MeCP2 shRNA')),
]

s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch')

# test_results = add_stat_annotation(ax, data=figure_data, x='condition', y='result',hue='group',
#                                    box_pairs=box_pairs,
#                                    text_annot_custom=["Glial marker", "Neuronal marker"],
#                                    perform_stat_test=False, pvalues=[0, 0],
#                                    loc='outside', verbose=0)

plt.xlabel("")
plt.ylabel("Relative Luciferase Activity")
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_3C.pdf', dpi=600, bbox_inches='tight')
#%%
# 3D is a gel band plot
# 3E
figure_data = data[data["figure"]=="3E"]

fig, ax = plt.subplots(figsize=(4,5))
ax = sns.barplot(x='cell_type', y='result', data=figure_data, hue="group",
                 ci="sd",
                 errwidth=1
                 )

# Add jitter with the swarmplot function
ax = sns.stripplot(x='cell_type', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none'
                   )
handles, labels = ax.get_legend_handles_labels()
# l = plt.legend(handles[-2:], labels[-2:], borderaxespad=0.)

box_pairs = [
    (('Neurons', 'Con shRNA'), ('Neurons', 'MeCP2 shRNA')),
    (('Astrocytes', 'Con shRNA'), ('Astrocytes', 'MeCP2 shRNA')),
]

s = add_stat_annotation(ax, data=figure_data, x='cell_type', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch')


plt.xlabel("")
plt.ylabel("Relative BDNF mRNA Expression")
ax.get_legend().remove()
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_3E.pdf', dpi=600, bbox_inches='tight')