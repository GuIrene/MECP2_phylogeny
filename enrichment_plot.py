import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.io as pio
from matplotlib.lines import Line2D
pio.renderers.default = "browser"
#%%
N_TO_SHOW = 10
tissues = pd.read_csv("/Users/ireneu/berman_lab/Rett/GeneAnalytics_tissue_cell.csv")
organs = pd.read_csv("/Users/ireneu/berman_lab/Rett/organsystem_results.txt", sep="\t", header=None, names=["Tissue / System", "Score"])
#%%
top_organs = organs.iloc[:N_TO_SHOW,:]
top_organ_names = top_organs["Tissue / System"].values
top_tissues = tissues.loc[[x in top_organ_names for x in tissues["Organ / Tissue"].values]]

#%%

fig, ax = plt.subplots(figsize=(6,4))
# Draw total score
ax = sns.barplot(data=top_organs, x="Tissue / System", y="Score",
                 color=(0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
                 ci=None,)

# Add jitter
ax = sns.stripplot(x='Organ / Tissue', y='Score', data=top_tissues, dodge=True,color='black', alpha=0.7,
                   jitter=True, marker='o', facecolors='none', order=top_organs["Tissue / System"].values,
                   )

for coll in ax.collections:
    y = coll.get_offsets()[0][1]
    coll.set_sizes([15*y])
ax.set_xticklabels(labels = top_organs["Tissue / System"].values, rotation=30)
plt.xticks(rotation=45)
plt.ylabel("System Score")

legend_elements = [Line2D([0], [0], marker='o', color='w', label=str(y),
                          markerfacecolor='black', markersize=1.5*y) for y in [2, 4, 6, 8]]

# Create the figure
ax.legend(title="Tissue Score", handles=legend_elements)

plt.show()
# plt.savefig('/Users/ireneu/berman_lab/Rett/jitter_plots/geneanalytics_390.pdf', dpi=600, bbox_inches='tight')
#%%
