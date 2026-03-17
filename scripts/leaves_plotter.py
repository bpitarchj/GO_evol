#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 12:04:53 2024

@author: bpitarch
"""
import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data_folder = os.listdir("data/")
data_folder.sort()
print(data_folder)


data =  [] #number of terms, new_terms, old_terms, obsolete_terms, new_obs, old_obs
years = []
for i in range(0,len(data_folder)):

    folder = data_folder[i]
    path = "data/"+folder+"/"    
    if os.path.isdir(path):
        years.append(int(folder))
        if i == 0:
            df = pd.read_csv(path+folder+"_leaf_terms",sep = "\t", header= None, names =["GO_term","Ontology_subtype","depth"])
            print(df)
            n_leaf= df.shape[0]
            n_leaf_BP = df[df["Ontology_subtype"]=="biological_process"].shape[0]
            n_leaf_MF = df[df["Ontology_subtype"]=="molecular_function"].shape[0]
            n_leaf_CC = df[df["Ontology_subtype"]=="cellular_component"].shape[0]
            
            max_depth = df["depth"].max()
            max_depth_BP = df[df["Ontology_subtype"]=="biological_process"]["depth"].max()
            max_depth_MF = df[df["Ontology_subtype"]=="molecular_function"]["depth"].max()
            max_depth_CC = df[df["Ontology_subtype"]=="cellular_component"]["depth"].max()
            
            mean_depth = df["depth"].mean()
            mean_depth_BP = df[df["Ontology_subtype"]=="biological_process"]["depth"].mean()
            mean_depth_MF = df[df["Ontology_subtype"]=="molecular_function"]["depth"].mean()
            mean_depth_CC = df[df["Ontology_subtype"]=="cellular_component"]["depth"].mean()
            
            new_leaf = np.nan
            new_leaf_BP = np.nan
            new_leaf_MF = np.nan
            new_leaf_CC = np.nan
            
            old_leaf = np.nan
            old_leaf_BP = np.nan
            old_leaf_MF = np.nan
            old_leaf_CC = np.nan   
            
            ###PLOT details
            ax = plt.figure(figsize = (35,10))
            plt.title("Distribution of depth in leaf terms in GO " + folder)
            ax = sns.kdeplot(df["depth"], color = "grey",label = "total")
            ax = sns.kdeplot(df[df["Ontology_subtype"]=="biological_process"]["depth"], color = "blue",label = "GO_BP")
            ax = sns.kdeplot(df[df["Ontology_subtype"]=="molecular_function"]["depth"], color = "black",label = "GO_MF")
            ax = sns.kdeplot(df[df["Ontology_subtype"]=="cellular_component"]["depth"], color = "green",label = "GO_CC")
            plt.legend()
            plt.savefig("results/"+folder+"_depth_distribution.jpg")
            
            depth_df_BP = df[df["Ontology_subtype"]=="biological_process"]["depth"].to_frame(name=folder)
            depth_df_MF = df[df["Ontology_subtype"]=="molecular_function"]["depth"].to_frame(name=folder)
            depth_df_CC = df[df["Ontology_subtype"]=="cellular_component"]["depth"].to_frame(name=folder)
            ###
            data.append([n_leaf,n_leaf_BP,n_leaf_MF,n_leaf_CC,
                         max_depth,max_depth_BP,max_depth_MF,max_depth_CC,
                         mean_depth,mean_depth_BP,mean_depth_MF,mean_depth_CC,
                         new_leaf, new_leaf_BP,new_leaf_MF,new_leaf_CC,
                         old_leaf, old_leaf_BP, old_leaf_MF,old_leaf_CC])
            print(data)
        else:
            df1 = pd.read_csv(path+folder+"_leaf_terms",sep = "\t", header= None, names =["GO_term","Ontology_subtype","depth"])
            print(df1)
            n_leaf1= df1.shape[0]
            n_leaf_BP1 = df1[df1["Ontology_subtype"]=="biological_process"].shape[0]
            n_leaf_MF1 = df1[df1["Ontology_subtype"]=="molecular_function"].shape[0]
            n_leaf_CC1 = df1[df1["Ontology_subtype"]=="cellular_component"].shape[0]
            
            max_depth1 = df1["depth"].max()
            max_depth_BP1 = df1[df1["Ontology_subtype"]=="biological_process"]["depth"].max()
            max_depth_MF1 = df1[df1["Ontology_subtype"]=="molecular_function"]["depth"].max()
            max_depth_CC1 = df1[df1["Ontology_subtype"]=="cellular_component"]["depth"].max()
            
            mean_depth1 = df1["depth"].mean()
            mean_depth_BP1 = df1[df1["Ontology_subtype"]=="biological_process"]["depth"].mean()
            mean_depth_MF1 = df1[df1["Ontology_subtype"]=="molecular_function"]["depth"].mean()
            mean_depth_CC1 = df1[df1["Ontology_subtype"]=="cellular_component"]["depth"].mean()
            
            new_leaf = len(set(list(df1["GO_term"]))-set(list(df["GO_term"])))
            new_leaf_BP = len(set(list(df1[df1["Ontology_subtype"]=="biological_process"]["GO_term"]))-set(list(df[df["Ontology_subtype"]=="biological_process"]["GO_term"])))
            new_leaf_MF = len(set(list(df1[df1["Ontology_subtype"]=="molecular_function"]["GO_term"]))-set(list(df[df["Ontology_subtype"]=="molecular_function"]["GO_term"])))
            new_leaf_CC = len(set(list(df1[df1["Ontology_subtype"]=="cellular_component"]["GO_term"]))-set(list(df[df["Ontology_subtype"]=="cellular_component"]["GO_term"])))
            
            old_leaf = len(set(list(df["GO_term"]))-set(list(df1["GO_term"])))
            old_leaf_BP = len(set(list(df[df["Ontology_subtype"]=="biological_process"]["GO_term"]))-set(list(df1[df1["Ontology_subtype"]=="biological_process"]["GO_term"])))
            old_leaf_MF = len(set(list(df[df["Ontology_subtype"]=="molecular_function"]["GO_term"]))-set(list(df1[df1["Ontology_subtype"]=="molecular_function"]["GO_term"])))
            old_leaf_CC = len(set(list(df[df["Ontology_subtype"]=="cellular_component"]["GO_term"]))-set(list(df1[df1["Ontology_subtype"]=="cellular_component"]["GO_term"])))
        
            
            ###PLOT details
            ax = plt.figure(figsize = (35,10))
            plt.title("Distribution of depth in leaf terms in GO " + folder)
            ax = sns.kdeplot(df1["depth"], color = "grey",label = "total")
            ax = sns.kdeplot(df1[df1["Ontology_subtype"]=="biological_process"]["depth"], color = "blue",label = "GO_BP")
            ax = sns.kdeplot(df1[df1["Ontology_subtype"]=="molecular_function"]["depth"], color = "black",label = "GO_MF")
            ax = sns.kdeplot(df1[df1["Ontology_subtype"]=="cellular_component"]["depth"], color = "green",label = "GO_CC")
            plt.legend()
            plt.savefig("results/"+folder+"_depth_distribution.jpg")
            
            
            depth_df_BP_1 = df[df["Ontology_subtype"]=="biological_process"]["depth"].to_frame(name=folder)
            depth_df_MF_1 = df[df["Ontology_subtype"]=="molecular_function"]["depth"].to_frame(name=folder)
            depth_df_CC_1 = df[df["Ontology_subtype"]=="cellular_component"]["depth"].to_frame(name=folder)
            depth_df_BP = pd.concat([depth_df_BP, depth_df_BP_1], axis = 1)
            depth_df_MF = pd.concat([depth_df_MF, depth_df_MF_1], axis = 1)
            depth_df_CC = pd.concat([depth_df_CC, depth_df_CC_1], axis = 1)
            ###
            
            data.append([n_leaf1,n_leaf_BP1,n_leaf_MF1,n_leaf_CC1,
                         max_depth1,max_depth_BP1,max_depth_MF1,max_depth_CC1,
                         mean_depth1,mean_depth_BP1,mean_depth_MF1,mean_depth_CC1,
                         new_leaf, new_leaf_BP,new_leaf_MF,new_leaf_CC,
                         old_leaf, old_leaf_BP, old_leaf_MF,old_leaf_CC])
            
            df = df1
columns = ["total","BP","MF","CC",",max_depth","max_depth_BP","max_depth_MF", "max_depth_CC",
           "mean_depth","mean_depth_BP","mean_depth_MF","mean_depth_CC",
           "new","new_BP","new_MF","new_CC","old","old_BP","old_MF","old_CC"]
data = pd.DataFrame(data, index= years, columns = columns)
print(data)

#%%
def Z_score(serie):
    return((serie-serie.mean())/serie.std())

def vline_plotter(serie):
    for value in list(serie):
        plt.axvline(value,color="grey", linewidth=1)
def leaf_norm(leaves,nodes):
    return((nodes-leaves)/leaves)

nodes_info = pd.read_csv("results/GO_evolution_data_nodes.csv", header = 0, index_col=0)
nodes= nodes_info[["terms_BP","terms_MF","terms_CC"]]
nodes.columns=["BP","MF","CC"]


#%%

#Plot evolution of terms
ax = plt.figure(figsize=(12,5))
plt.rc("font",size=16)
#plt.title("Ratio of no leaf to leaf terms in GO")
plt.xlabel("Year")
plt.ylabel("No leaf terms / Leaf terms")
ax=plt.plot(data.index,leaf_norm(data["BP"],nodes["BP"]),label="GO:BP",marker="o" ,linewidth =  2.5 , color = "blue")
ax=plt.plot(data.index,leaf_norm(data["MF"],nodes["MF"]),label="GO:MF",marker="o" , linewidth =  2.5, color="black")
ax=plt.plot(data.index,leaf_norm(data["CC"],nodes["CC"]),label="GO:CC", marker="o", linewidth =  2.5, color="green")
plt.xticks(years, rotation=315)
vline_plotter(data.index)
plt.legend(loc="upper left")
plt.tight_layout()
#plt.savefig("results/GO_evolution_number_of_leaf_terms.jpg")
plt.savefig("new_figures_monica/GO_evolution_number_of_leaf_terms.svg")

#%%
#Plot evolution of depth

ax = plt.figure(figsize = (15,8))
#plt.title("Depth of GO:BP terms")
sns.violinplot(data=depth_df_BP , color="blue",linewidth=1)
plt.rc("font",size = 16)
plt.ylabel("Depth of terms")
plt.xlabel("Year")
plt.tick_params("x",rotation = 315)
plt.tight_layout()
#plt.savefig("results/GO-BP_evolution_depth_terms.jpg",dpi = 200)
plt.savefig("new_figures_monica/GO-BP_evolution_depth_terms.svg",dpi = 200)
plt.close()

ax = plt.figure(figsize = (15,8))
plt.title("Depth of GO:MF terms")
sns.violinplot(data=depth_df_MF , color="dimgrey",linewidth=1)
plt.rc("font",size = 16)
plt.ylabel("Depth of terms")
plt.xlabel("Year")
plt.tick_params("x",rotation = 315)
plt.tight_layout()
plt.savefig("results/GO-MF_evolution_depth_terms.jpg",dpi = 200)
plt.close()

ax = plt.figure(figsize = (15,8))
plt.title("Depth of GO:CC terms")
sns.violinplot(data=depth_df_CC , color="green",linewidth=1)
plt.rc("font",size = 16)
plt.ylabel("Depth of terms")
plt.xlabel("Year")
plt.tick_params("x",rotation = 315)
plt.tight_layout()
plt.savefig("results/GO-CC_evolution_depth_terms.jpg",dpi = 200)
plt.close()




#%%
"""
#Plot evolution of new terms

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of number of new leaf GO terms",size = 100)
plt.xlabel("Year")
axs[0].title.set_text("Number of new leaf GO terms")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of new leaf terms")
axs[1].title.set_text("Number of new leaf GO terms (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["new_BP"]-data["new_BP"].mean())/data["new_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["new_MF"]-data["new_MF"].mean())/data["new_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["new_CC"]-data["new_CC"].mean())/data["new_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of new leaf terms")

plt.savefig("results/GO_evolution_number_of_new_leaf_terms.svg")
#%%

#Plot evolution of old terms

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of number of old leaf GO terms",size = 100)
plt.xlabel("Year")
axs[0].title.set_text("Number of old leaf GO terms")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of old leaf terms")
axs[1].title.set_text("Number of old leaf GO terms (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["old_BP"]-data["old_BP"].mean())/data["old_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["old_MF"]-data["old_MF"].mean())/data["old_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["old_CC"]-data["old_CC"].mean())/data["old_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of old leaf terms")

plt.savefig("results/GO_evolution_number_of_old_leaf_terms.svg")

#%%

#Plot evolution of depth

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of max depth leaf GO terms",size = 100)
plt.xlabel("Year")
axs[0].title.set_text("Max depth of leaf GO terms")
axs[0] = sns.lineplot(data=data, x =data.index, y=data["max_depth_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5, color="blue")
axs[0] = sns.lineplot(data=data, x =data.index, y=data["max_depth_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0] = sns.lineplot(data=data, x =data.index, y=data["max_depth_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("Max depth leaf terms")

axs[1].title.set_text("Mean depth leaf GO terms")
axs[1] = sns.lineplot(data=data, x =data.index, y=data["mean_depth_BP"],label="GO_BP", ax=axs[1], linewidth =  2.5, color="blue")
axs[1] = sns.lineplot(data=data, x =data.index, y=data["mean_depth_MF"],label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1] = sns.lineplot(data=data, x =data.index, y=data["mean_depth_CC"],label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("Mean depth leaf terms")


plt.savefig("results/GO_evolution_depth_leaf_terms.svg")
"""
#%%

#Save data
data.to_csv("results/GO_evolution_leaf_terms.csv",header = True, index = True)
#%%

### Cross with historical series

