import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom

#%%
all_genes = pd.read_csv("/Users/ireneu/berman_lab/Rett/protein_coding_biomart.txt", sep = "\t")
chromosomes = [str(x) for x in range(1,23)]+["X", "Y", "MT"]
chrom_filter =  [x in chromosomes for x in all_genes["Chromosome/scaffold name"].values]
all_genes = all_genes.iloc[chrom_filter,:]

#%%
query = pd.read_csv("/Users/ireneu/berman_lab/Rett/390_eu_m.txt")
query.at[-1] = "MECP2"
#%%
background = pd.read_csv("/Users/ireneu/berman_lab/Rett/20192_background_set.txt")
#%%
merged = pd.merge(background, all_genes, how="left", left_on="Gene", right_on="Gene name")
merged.dropna(subset=["Chromosome/scaffold name"],inplace=True)
merged.drop_duplicates(subset=["Gene"], keep="first", inplace=True)
merged["co-evolved"] = merged["Gene"].apply(lambda x: x in query["Gene"].values)
#%%

population = len(merged)
sample_size = sum(merged["co-evolved"])

#%%
all_pvals = []
names = []
observed = []
expected = []
grouped = merged.groupby(["Chromosome/scaffold name", "Karyotype band"])
for group_name in grouped.groups.keys():
    if group_name[0] == "MT":
        continue
    df = grouped.get_group(group_name)
    names.append(group_name[0]+str(group_name[1]))
    total_in_loc = len(df)
    expected.append(total_in_loc/population*sample_size)
    positive_in_loc = sum(df["co-evolved"])
    observed.append(positive_in_loc)
    pval = hypergeom.sf(positive_in_loc - 1, population, total_in_loc, sample_size)
    all_pvals.append(pval)

#%%
sorting_order = np.argsort(observed)[::-1]

sns.set_style("white")
n_to_show = 20
fig, ax = plt.subplots()
x = np.arange(n_to_show)  # the label locations
width = 0.5  # the width of the bars

fig, ax = plt.subplots(figsize=(14,10))
rects1 = ax.bar(x - width/2, np.array(expected)[sorting_order][:n_to_show],
                width, label='expected')
rects2 = ax.bar(x + width/2, np.array(observed)[sorting_order][:n_to_show], width, label='observed/expected')
for i in x:
    if np.array(all_pvals)[sorting_order][i] < 0.001:
        ax.text(x=x[i] -0.1, y=np.array(observed)[sorting_order][i], s='***', fontdict={"size":20})
    elif np.array(all_pvals)[sorting_order][i] < 0.01:
        ax.text(x=x[i], y=np.array(observed)[sorting_order][i], s='**', fontdict={"size":20})
    elif np.array(all_pvals)[sorting_order][i] < 0.05:
        ax.text(x=x[i]+0.1, y=np.array(observed)[sorting_order][i], s='*', fontdict={"size":20})
ax.set_xticks(x)
ax.set_xticklabels(np.array(names)[sorting_order][:n_to_show], rotation=45, size=20)
plt.ylabel("Number of Genes", size=20)
plt.xlabel("Chromosomal bands (ranked by observed count)", size = 20)

# plt.savefig("/Users/ireneu/berman_lab/Rett/391_genes_karyotype_band_by_observed.pdf")
plt.show()