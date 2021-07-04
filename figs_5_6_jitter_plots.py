import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
#%%
input_path = "/Users/ireneu/berman_lab/Rett/Figure5_raw_results.csv"
data = pd.read_csv(input_path)
#%%
color_palette = [(0.12156862745098039, 0.4666666666666667, 0.7058823529411765),"grey",
 (1.0, 0.4980392156862745, 0.054901960784313725), "#fcc603"]
hue_order=['GFAP Con shRNA','bIII-tubulin Con shRNA','GFAP MeCP2 shRNA','bIII-tubulin MeCP2 shRNA']

#5A
figure_data = data[data["figure"]=="5A"]

fig, ax = plt.subplots(figsize=(8,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group", palette=color_palette,
                 ci="sd",
                 errwidth=1,
                 hue_order=hue_order
                 )

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none', hue_order=hue_order

                   )
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-4:], labels[-4:], loc=1, borderaxespad=0.)

box_pairs = [
    # (('Control', 'GFAP Con shRNA'), ('Control', 'GFAP MeCP2 shRNA')),
    # (('Control', 'bIII-tubulin Con shRNA'), ('Control', 'bIII-tubulin MeCP2 shRNA')),

    # (('Pacritinib', 'GFAP Con shRNA'),('Control', 'GFAP Con shRNA')),
    (('Pacritinib', 'GFAP MeCP2 shRNA'), ('Control', 'GFAP MeCP2 shRNA')),
    # (('Pacritinib', 'bIII-tubulin Con shRNA'), ('Control', 'bIII-tubulin Con shRNA')),
    (('Pacritinib', 'bIII-tubulin MeCP2 shRNA'), ('Control', 'bIII-tubulin MeCP2 shRNA')),

    # (('DMF', 'GFAP Con shRNA'), ('Control', 'GFAP Con shRNA')),
    (('DMF', 'GFAP MeCP2 shRNA'), ('Control', 'GFAP MeCP2 shRNA')),
    # (('DMF', 'bIII-tubulin Con shRNA'), ('Control', 'bIII-tubulin Con shRNA')),
    (('DMF', 'bIII-tubulin MeCP2 shRNA'), ('Control', 'bIII-tubulin MeCP2 shRNA')),

    # (('EPO', 'GFAP Con shRNA'), ('Control', 'GFAP Con shRNA')),
    (('EPO', 'GFAP MeCP2 shRNA'), ('Control', 'GFAP MeCP2 shRNA')),
    # (('EPO', 'bIII-tubulin Con shRNA'), ('Control', 'bIII-tubulin Con shRNA')),
    (('EPO', 'bIII-tubulin MeCP2 shRNA'), ('Control', 'bIII-tubulin MeCP2 shRNA')),
]

s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch',
                        hue_order=hue_order
                        )



plt.xlabel("")
plt.ylabel("Relative mRNA Expression")
ax.get_legend().remove()
# plt.show()

plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_5A.pdf', dpi=600, bbox_inches='tight')
#%%
#5B
figure_data = data[data["figure"]=="5B"]

fig, ax = plt.subplots(figsize=(5,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group",
                 ci="sd",
                 errwidth=1
                 )

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none'
                   )
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-2:], labels[-2:], loc=4, borderaxespad=0.)

box_pairs = [
    # (('Control', 'Con shRNA'), ('Control', 'MeCP2 shRNA')),
    #
    # (('Control', 'Con shRNA'), ('Pacritinib', 'Con shRNA')),
    # (('Control', 'Con shRNA'), ('DMF', 'Con shRNA')),
    # (('Control', 'Con shRNA'), ('EPO', 'Con shRNA')),
    #
    (('Control', 'MeCP2 shRNA'), ('Pacritinib', 'MeCP2 shRNA')),
    (('Control', 'MeCP2 shRNA'), ('DMF', 'MeCP2 shRNA')),
    (('Control', 'MeCP2 shRNA'), ('EPO', 'MeCP2 shRNA')),

]

s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch')


plt.xlabel("")
plt.ylabel("Relative EAAT2 Expression")
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_5B.pdf', dpi=600, bbox_inches='tight')
#%%
input_path = "/Users/ireneu/berman_lab/Rett/Figure6_raw_results.csv"
data = pd.read_csv(input_path)
#6A
figure_data = data[data["figure"]=="6A"]

fig, ax = plt.subplots(figsize=(5,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group",
                 ci="sd",
                 errwidth=1
                 )


# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none'
                   )
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-2:], labels[-2:], loc='lower right', borderaxespad=0.)

box_pairs = [
    # (('Control', 'Con shRNA'), ('Control', 'MeCP2 shRNA')),

    # (('Pacritinib', 'Con shRNA'), ('Control', 'Con shRNA')),
    # (('DMF', 'Con shRNA'), ('Control', 'Con shRNA')),
    # (('EPO', 'Con shRNA'), ('Control', 'Con shRNA')),

    (('Pacritinib', 'MeCP2 shRNA'), ('Control', 'MeCP2 shRNA')),
    (('DMF', 'MeCP2 shRNA'), ('Control', 'MeCP2 shRNA')),
    (('EPO', 'MeCP2 shRNA'), ('Control', 'MeCP2 shRNA')),

]

s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch')


plt.xlabel("")
plt.ylabel("Relative BDNF mRNA Expression")
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_6A.pdf', dpi=600, bbox_inches='tight')
#%%
#6B
figure_data = data[data["figure"]=="6B"]

fig, ax = plt.subplots(figsize=(5,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group",
                 ci="sd",
                 errwidth=1
                 )

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none'
                   )
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-2:], labels[-2:], loc=4, borderaxespad=0.)

box_pairs = [
    # (('Control', 'Con shRNA'), ('Control', 'MeCP2 shRNA')),

    # (('Control', 'Con shRNA'), ('Pacritinib', 'Con shRNA')),
    # (('Control', 'Con shRNA'), ('DMF', 'Con shRNA')),
    # (('Control', 'Con shRNA'), ('EPO', 'Con shRNA')),

    (('Control', 'MeCP2 shRNA'), ('Pacritinib', 'MeCP2 shRNA')),
    (('Control', 'MeCP2 shRNA'), ('DMF', 'MeCP2 shRNA')),
    (('Control', 'MeCP2 shRNA'), ('EPO', 'MeCP2 shRNA')),

]

s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch')


plt.xlabel("")
plt.ylabel("BDNF Secretion - pg/ml")
plt.ylim((0,50))

# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_6B.pdf', dpi=600, bbox_inches='tight')
#%%
#6C
figure_data = data[data["figure"]=="6C"]

fig, ax = plt.subplots(figsize=(5,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group",
                 ci="sd",
                 errwidth=1
                 )

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none'
                   )
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-2:], labels[-2:], loc=4, borderaxespad=0.)

box_pairs = [
    # (('Control', 'Con shRNA'), ('Control', 'MeCP2 shRNA')),

    # (('Control', 'Con shRNA'), ('DMF+EPO', 'Con shRNA')),
    # (('Control', 'Con shRNA'), ('DMF', 'Con shRNA')),
    # (('Control', 'Con shRNA'), ('EPO', 'Con shRNA')),

    (('Control', 'MeCP2 shRNA'), ('DMF+EPO', 'MeCP2 shRNA')),
    # (('Control', 'MeCP2 shRNA'), ('DMF', 'MeCP2 shRNA')),
    # (('Control', 'MeCP2 shRNA'), ('EPO', 'MeCP2 shRNA')),

    (('DMF+EPO', 'MeCP2 shRNA'), ('DMF', 'MeCP2 shRNA')),
    (('DMF+EPO', 'MeCP2 shRNA'), ('EPO', 'MeCP2 shRNA'))

]


s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch')


plt.xlabel("")
plt.ylim((0,50))
plt.ylabel("BDNF Secretion - pg/ml")
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_6C.pdf', dpi=600, bbox_inches='tight')