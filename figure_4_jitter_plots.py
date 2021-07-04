import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
#%%
input_path = "/Users/ireneu/berman_lab/Rett/Figure4_raw_results.csv"
data = pd.read_csv(input_path)
#%%
color_palette = [(0.12156862745098039, 0.4666666666666667, 0.7058823529411765),"grey",
 (1.0, 0.4980392156862745, 0.054901960784313725), "#fcc603"]
hue_order=['PBS Control shRNA','EPO Control shRNA','PBS MeCP2 shRNA','EPO MeCP2 shRNA']
#4A
figure_data = data[data["figure"]=="4A"]

fig, ax = plt.subplots(figsize=(8,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group", palette=color_palette, ci='sd',
                 errwidth=1, hue_order=hue_order)

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color='none', hue="group", dodge=True,
                   facecolors='none', jitter=True, hue_order=hue_order)
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-4:], labels[-4:], loc=2, borderaxespad=0.)

box_pairs = [
    # (('CD206', 'PBS Control shRNA'), ('CD206', 'PBS MeCP2 shRNA')),
    # (('CD206', 'PBS Control shRNA'), ('CD206', 'EPO Control shRNA')),
    (('CD206', 'EPO MeCP2 shRNA'), ('CD206', 'PBS MeCP2 shRNA')),

    # (('IL-13', 'PBS Control shRNA'), ('IL-13', 'PBS MeCP2 shRNA')),
    # (('IL-13', 'PBS Control shRNA'), ('IL-13', 'EPO Control shRNA')),
    (('IL-13','EPO MeCP2 shRNA'), ('IL-13', 'PBS MeCP2 shRNA')),

    # (('CD86', 'PBS Control shRNA'), ('CD86', 'PBS MeCP2 shRNA')),
    # (('CD86', 'PBS Control shRNA'), ('CD86', 'EPO Control shRNA')),
    (('CD86', 'EPO MeCP2 shRNA'), ('CD86', 'PBS MeCP2 shRNA')),

    # (('IL-1', 'PBS Control shRNA'), ('IL-1', 'PBS MeCP2 shRNA')),
    # (('IL-1', 'PBS Control shRNA'), ('IL-1', 'EPO Control shRNA')),
    (('IL-1', 'EPO MeCP2 shRNA'), ('IL-1', 'PBS MeCP2 shRNA')),

    # (('TNF-alpha', 'PBS Control shRNA'), ('TNF-alpha', 'PBS MeCP2 shRNA')),
    # (('TNF-alpha', 'PBS Control shRNA'), ('TNF-alpha', 'EPO Control shRNA')),
    (('TNF-alpha', 'EPO MeCP2 shRNA'), ('TNF-alpha', 'PBS MeCP2 shRNA')),

]
s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch',
                        hue_order=hue_order)

# test_results = add_stat_annotation(ax, data=figure_data, x='condition', y='result',
#                                    box_pairs=[('CD206', 'IL-13'),('CD86', 'TNF-alpha')],
#                                    text_annot_custom=["M2 markers", "M1 markers"],
#                                    perform_stat_test=False, pvalues=[0, 0],
#                                    loc='outside', verbose=0)

plt.xlabel("")
plt.ylabel("Relative mRNA Expression")
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_4A.pdf', dpi=600, bbox_inches='tight')
#%%
hue_order=['PBS Control shRNA','DMF Control shRNA','PBS MeCP2 shRNA','DMF MeCP2 shRNA']


#4B
figure_data = data[data["figure"]=="4B"]

fig, ax = plt.subplots(figsize=(8,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group", palette=color_palette,
                 ci='sd',
                 errwidth=1,
                 hue_order=hue_order
                 )

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True,marker='o',facecolors='none', hue_order=hue_order)
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-4:], labels[-4:], loc=2, borderaxespad=0.)

box_pairs = [
    # (('CD206', 'PBS Control shRNA'), ('CD206', 'PBS MeCP2 shRNA')),
    # (('CD206', 'PBS Control shRNA'), ('CD206', 'DMF Control shRNA')),
    (('CD206', 'DMF MeCP2 shRNA'), ('CD206', 'PBS MeCP2 shRNA')),

    # (('IL-13', 'PBS Control shRNA'), ('IL-13', 'PBS MeCP2 shRNA')),
    # (('IL-13', 'PBS Control shRNA'), ('IL-13', 'DMF Control shRNA')),
    (('IL-13','DMF MeCP2 shRNA'), ('IL-13', 'PBS MeCP2 shRNA')),

    # (('CD86', 'PBS Control shRNA'), ('CD86', 'PBS MeCP2 shRNA')),
    # (('CD86', 'PBS Control shRNA'), ('CD86', 'DMF Control shRNA')),
    (('CD86', 'DMF MeCP2 shRNA'), ('CD86', 'PBS MeCP2 shRNA')),

    # (('IL-1', 'PBS Control shRNA'), ('IL-1', 'PBS MeCP2 shRNA')),
    # (('IL-1', 'PBS Control shRNA'), ('IL-1', 'DMF Control shRNA')),
    (('IL-1', 'DMF MeCP2 shRNA'), ('IL-1', 'PBS MeCP2 shRNA')),

    # (('TNF-alpha', 'PBS Control shRNA'), ('TNF-alpha', 'PBS MeCP2 shRNA')),
    # (('TNF-alpha', 'PBS Control shRNA'), ('TNF-alpha', 'DMF Control shRNA')),
    (('TNF-alpha', 'DMF MeCP2 shRNA'), ('TNF-alpha', 'PBS MeCP2 shRNA')),

]
s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch',
                        hue_order=hue_order
                        )

# test_results = add_stat_annotation(ax, data=figure_data, x='condition', y='result',
#                                    box_pairs=[('CD206', 'IL-13'),('CD86', 'TNF-alpha')],
#                                    text_annot_custom=["M2 markers", "M1 markers"],
#                                    perform_stat_test=False, pvalues=[0, 0],
#                                    loc='outside', verbose=0)

plt.xlabel("")
plt.ylabel("Relative mRNA Expression")
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_4B.pdf', dpi=600, bbox_inches='tight')
#%%
#4C
figure_data = data[data["figure"]=="4C"]
hue_order=['PBS Control shRNA','Pacritinib Control shRNA','PBS MeCP2 shRNA','Pacritinib MeCP2 shRNA']

fig, ax = plt.subplots(figsize=(8,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group", palette=color_palette,
                 ci='sd',
                 errwidth=1,
                 hue_order=hue_order

                 )

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none', hue_order=hue_order

                   )
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-4:], labels[-4:], loc=2, borderaxespad=0.)


box_pairs = [
    # (('CD206', 'PBS Control shRNA'), ('CD206', 'PBS MeCP2 shRNA')),
    # (('CD206', 'PBS Control shRNA'), ('CD206', 'Pacritinib Control shRNA')),
    (('CD206', 'Pacritinib MeCP2 shRNA'), ('CD206', 'PBS MeCP2 shRNA')),

    # (('IL-13', 'PBS Control shRNA'), ('IL-13', 'PBS MeCP2 shRNA')),
    # (('IL-13', 'PBS Control shRNA'), ('IL-13', 'Pacritinib Control shRNA')),
    (('IL-13','Pacritinib MeCP2 shRNA'), ('IL-13', 'PBS MeCP2 shRNA')),

    # (('CD86', 'PBS Control shRNA'), ('CD86', 'PBS MeCP2 shRNA')),
    # (('CD86', 'PBS Control shRNA'), ('CD86', 'Pacritinib Control shRNA')),
    (('CD86', 'Pacritinib MeCP2 shRNA'), ('CD86', 'PBS MeCP2 shRNA')),

    # (('IL-1', 'PBS Control shRNA'), ('IL-1', 'PBS MeCP2 shRNA')),
    # (('IL-1', 'PBS Control shRNA'), ('IL-1', 'Pacritinib Control shRNA')),
    (('IL-1', 'Pacritinib MeCP2 shRNA'), ('IL-1', 'PBS MeCP2 shRNA')),

    # (('TNF-alpha', 'PBS Control shRNA'), ('TNF-alpha', 'PBS MeCP2 shRNA')),
    # (('TNF-alpha', 'PBS Control shRNA'), ('TNF-alpha', 'Pacritinib Control shRNA')),
    (('TNF-alpha', 'Pacritinib MeCP2 shRNA'), ('TNF-alpha', 'PBS MeCP2 shRNA')),

]
s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch',
                        hue_order=hue_order
                        )

# test_results = add_stat_annotation(ax, data=figure_data, x='condition', y='result',
#                                    box_pairs=[('CD206', 'IL-13'),('CD86', 'TNF-alpha')],
#                                    text_annot_custom=["M2 markers", "M1 markers"],
#                                    perform_stat_test=False, pvalues=[0, 0],
#                                    loc='outside', verbose=0)

plt.xlabel("")
plt.ylabel("Relative mRNA Expression")
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_4C.pdf', dpi=600, bbox_inches='tight')
#%%
#4D
figure_data = data[data["figure"]=="4D"]

fig, ax = plt.subplots(figsize=(5,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group",
                 ci='sd',
                 errwidth=1
                 )

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none'
                   )
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-2:], labels[-2:], loc=4, borderaxespad=0.)

box_pairs = [
    # (('PBS', 'Con shRNA'), ('PBS', 'MeCP2 shRNA')),

    # (('PBS', 'Con shRNA'), ('Pacritinib', 'Con shRNA')),
    # (('PBS', 'Con shRNA'), ('DMF', 'Con shRNA')),
    # (('PBS', 'Con shRNA'), ('EPO', 'Con shRNA')),

    (('PBS', 'MeCP2 shRNA'), ('Pacritinib', 'MeCP2 shRNA')),
    (('PBS', 'MeCP2 shRNA'), ('DMF', 'MeCP2 shRNA')),
    (('PBS', 'MeCP2 shRNA'), ('EPO', 'MeCP2 shRNA')),
]

s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch')


plt.xlabel("")
plt.ylabel("% Phagocytosis/Control")
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_4D.pdf', dpi=600, bbox_inches='tight')
#%%
#4E
figure_data = data[data["figure"]=="4E"]

fig, ax = plt.subplots(figsize=(5,5))
ax = sns.barplot(x='condition', y='result', data=figure_data, hue="group",
                 ci='sd',
                 errwidth=1
                 )

# Add jitter with the swarmplot function
ax = sns.stripplot(x='condition', y='result', data=figure_data, color="black", hue="group", dodge=True,
                   jitter=True, marker='o', facecolors='none'
                   )
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[-2:], labels[-2:], loc=4, borderaxespad=0.)

box_pairs = [
    # (('PBS', 'Con shRNA'), ('PBS', 'MeCP2 shRNA')),

    # (('PBS', 'Con shRNA'), ('Pacritinib', 'Con shRNA')),
    # (('PBS', 'Con shRNA'), ('DMF', 'Con shRNA')),
    # (('PBS', 'Con shRNA'), ('EPO', 'Con shRNA')),

    (('PBS', 'MeCP2 shRNA'), ('Pacritinib', 'MeCP2 shRNA')),
    (('PBS', 'MeCP2 shRNA'), ('DMF', 'MeCP2 shRNA')),
    (('PBS', 'MeCP2 shRNA'), ('EPO', 'MeCP2 shRNA')),
]

s = add_stat_annotation(ax, data=figure_data, x='condition', y='result', hue='group', loc = 'inside',
                        box_pairs=box_pairs, comparisons_correction=None, text_format='star', test='t-test_welch')


plt.xlabel("")
plt.ylabel("Relative NF-kB Dependent Luciferase Activity")
# plt.show()
plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/fig_4E.pdf', dpi=600, bbox_inches='tight')