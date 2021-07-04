import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import upsetplot

#%%
data_dir = "/Users/ireneu/berman_lab/Rett/revisions/string_data/"
#load data
with gzip.open(data_dir+"9606.protein.links.detailed.v11.0.txt.gz") as infile:
    links = pd.read_csv(infile, sep=" ")
with gzip.open(data_dir+"9606.protein.info.v11.0.txt.gz") as infile:
    identifiers = pd.read_csv(infile, sep="\t")
#%%
id_to_name = dict(zip(identifiers["protein_external_id"], identifiers["preferred_name"]))
name_to_id = dict(zip(identifiers["preferred_name"], identifiers["protein_external_id"]))

#Some gene are listed in STRING under a different alias
extra_name_to_id = {"CFAP157":"9606.ENSP00000362392",
                    "RYDEN":"9606.ENSP00000253110",
                    "TRMO":"9606.ENSP00000364260",
                    "MFSD4A":"9606.ENSP00000356115",
                    "SHTN1":"9606.ENSP00000347532",
                    "WASHC4":"9606.ENSP00000484713",
                    "SUSD6":"9606.ENSP00000344424",
                    "RCC1L":"9606.ENSP00000480364",
                    "BORCS5":"9606.ENSP00000321546",
                    "SPATA46":"9606.ENSP00000356912",
                    "INO80B-WBP1":1, #not in STRING
                    "MROH7-TTC4":2,
                    "STRCP1":3,
                    "PLPBP":4,
                    "SEC11B":5
                    }
name_to_id.update(extra_name_to_id)
#%%
mecp2_links = links[links["protein1"] == name_to_id["MECP2"]] #information is double
#%%
eukaryotes = pd.read_csv("/Users/ireneu/berman_lab/Rett/eukaryotes_mecp2_correlation.csv")
mammals = pd.read_csv("/Users/ireneu/berman_lab/Rett/mammals_mecp2_correlation.csv")

eukaryotes["protein_id"] = eukaryotes["Gene"].apply(lambda x: name_to_id[x] if x in name_to_id.keys() else None)
mammals["protein_id"] = mammals["Gene"].apply(lambda x: name_to_id[x] if x in name_to_id.keys() else None)

#%%
eukaryotes.columns = ["gene", "PP_eukaryotes", "protein2"]
mammals.columns = ["gene", "PP_mammals", "protein2"]

eukaryotes = eukaryotes.iloc[1:201,:]
mammals = mammals.iloc[1:201,:]
pp = pd.merge(mammals, eukaryotes, on=["gene", "protein2"], how="outer")
pp.fillna(0, inplace=True)
pp["PP"] = pp.apply(lambda x: max(x["PP_eukaryotes"], x["PP_mammals"]), axis=1)
mecp2_links = pd.merge(mecp2_links, pp[["gene", "protein2", "PP"]], on="protein2", how="outer")

mecp2_links["protein1"]= mecp2_links["protein1"].fillna('9606.ENSP00000395535')
mecp2_links.fillna(0, inplace=True)

#%%

mecp2_links = pd.merge(mecp2_links, pp, on="protein2", how="outer")
mecp2_links["protein1"]= mecp2_links["protein1"].fillna('9606.ENSP00000395535')
mecp2_links.fillna(0, inplace=True)


#%%
binary_links = mecp2_links.dropna().drop(["gene", "combined_score", "cooccurence", "fusion", "neighborhood"], axis=1)
binary_links.iloc[:,2:] = binary_links.iloc[:,2:] > 0
binary_links["sum"]=np.sum(binary_links.iloc[:,2:], axis=1)
binary_links = binary_links[binary_links["sum"]>0]
binary_links.drop("sum", inplace=True, axis=1)
#%%
cols = binary_links.columns.difference(['protein1', 'protein2']).tolist()
s = binary_links.groupby(cols).size()
upsetplot.plot(s, show_counts='%d', sort_by='cardinality')
plt.show()
# plt.savefig("/Users/ireneu/berman_lab/Rett/390_genes_string_upsetplot.pdf")