import numpy as np
import pandas as pd
from IPython.display import display
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import seaborn as sns
from scipy.stats import pearsonr
import json
from itertools import combinations
import os
from scipy.stats import zscore
from matplotlib_venn import venn2

#%%

def correlations_to_gene(npp, gene):
    '''
    npp: Phylogenetic profile dataframe, with genes
    as index and species as columns
    gene: query
    '''
    res = pd.DataFrame()
    res["Gene"] = npp.index
    res["Correlation"] = res["Gene"].apply(lambda x:
             pearsonr(npp.loc[gene].values, npp.loc[x].values)[0])
    res.sort_values(by = "Correlation", ascending = False,
                           inplace = True)
    res.reset_index(inplace = True, drop = True)
    return res


#%%
NPP_PATH = "9606_NPP.tsv"
PHYLO_FILES_PATH = "phylo_files"
N_MAMMALS = 51
npp = pd.read_csv(NPP_PATH, sep = "\t") #load phylogenetic matrix
npp.set_index("Unnamed: 0", inplace = True)
mammalian_npp = npp.iloc[:,:N_MAMMALS]

#%%
all_eukaryotes = correlations_to_gene(npp, "MECP2")
mammalians = correlations_to_gene(mammalian_npp, "MECP2")
#%%
#Draw ranked correlations for all human genes
fig, ax = plt.subplots()
ax.step(np.arange(len(all_eukaryotes)), all_eukaryotes["Correlation"],
        label= "Eukaryotes")
ax.step(np.arange(len(mammalians)), mammalians["Correlation"],
       label = "Mammals", c="#E27C60")
ax.axvline(x = 200, c="black")
ax.set_title("Gene NPP Correlations with MECP2")
ax.set_ylabel("Pearson Correlation")
ax.set_xlabel("Human Genes")
ax.legend()
plt.show()
#%%
#Draw M200 and E200 intersection, fig1 C
top_200_m = mammalians["Gene"].iloc[:201].values
top_200_eu = all_eukaryotes["Gene"].iloc[:201].values

eu = set(top_200_eu)
m = set(top_200_m)
eu.remove("MECP2")
m.remove("MECP2")
# Second way
venn2([eu, m], set_colors=["#55BCC9", "#E27D60"],
      set_labels=["Eukaryotes", "Mammals"])
plt.show()
#%%
