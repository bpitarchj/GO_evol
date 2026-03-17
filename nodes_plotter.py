#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 15:47:44 2024

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
        years.append(folder)
        print(folder)
        if i == 0:
            df = pd.read_csv(path+folder+"_nodes",sep = "\t", header= None, names =["GO_term","Ontology_subtype"])
            print(df)
            nodes = df.shape[0]
            nodes_BP = df[df["Ontology_subtype"]=="biological_process"].shape[0]
            nodes_MF = df[df["Ontology_subtype"]=="molecular_function"].shape[0]
            nodes_CC = df[df["Ontology_subtype"]=="cellular_component"].shape[0]
            print(nodes_BP,nodes_MF,nodes_CC)
            new_nodes = np.nan
            new_nodes_BP = np.nan
            new_nodes_MF = np.nan
            new_nodes_CC = np.nan
            old_nodes = np.nan
            old_nodes_BP = np.nan
            old_nodes_MF = np.nan
            old_nodes_CC = np.nan
            obs_df = pd.read_csv(path+folder+"_obsolete_nodes",sep = "\t", header= None, names =["GO_term","Ontology_subtype"])
            obs_nodes = obs_df.shape[0]
            obs_nodes_BP = obs_df[obs_df["Ontology_subtype"]=="biological_process"].shape[0]
            obs_nodes_MF = obs_df[obs_df["Ontology_subtype"]=="molecular_function"].shape[0]
            obs_nodes_CC = obs_df[obs_df["Ontology_subtype"]=="cellular_component"].shape[0]
            new_obs_nodes = np.nan
            new_obs_nodes_BP = np.nan
            new_obs_nodes_MF = np.nan
            new_obs_nodes_CC = np.nan
            old_obs_nodes = np.nan
            old_obs_nodes_BP = np.nan
            old_obs_nodes_MF = np.nan
            old_obs_nodes_CC = np.nan
            
            historical_series = df
            historical_series["start"] = int(folder)
            historical_series["end"] = 2025
            print(historical_series)
            inter = list(set(list(df["GO_term"])) & set(list(obs_df["GO_term"])))
            if len(inter) != 0:
                print(folder)
                sys.exit()
            for term in list(obs_df["GO_term"]):
                row = historical_series[historical_series["GO_term"] == term]
                if row.shape[0] == 0:
                    continue
                elif row.shape[0] > 1:
                    print(row)
                    print("Error")
                    print(sys.exit())
                else:
                    index = row.index
                    historical_series.iloc[index,3] = int(folder)
            
            print(historical_series)
            data.append([nodes,nodes_BP,nodes_MF,nodes_CC,
                        new_nodes,new_nodes_BP,new_nodes_MF,new_nodes_CC,
                        old_nodes,old_nodes_BP,old_nodes_MF,old_nodes_CC,
                        obs_nodes,obs_nodes_BP,obs_nodes_MF,obs_nodes_CC,
                        new_obs_nodes,new_obs_nodes_BP,new_obs_nodes_MF,new_obs_nodes_CC,
                        old_obs_nodes,old_obs_nodes_BP,old_obs_nodes_MF,old_obs_nodes_CC])
            print(data)
        else:
            df1 = pd.read_csv(path+folder+"_nodes",sep = "\t", header= None, names =["GO_term","Ontology_subtype"])
            nodes1 = df1.shape[0]
            nodes_BP_1 = df1[df1["Ontology_subtype"]=="biological_process"].shape[0]
            nodes_MF_1 = df1[df1["Ontology_subtype"]=="molecular_function"].shape[0]
            nodes_CC_1 = df1[df1["Ontology_subtype"]=="cellular_component"].shape[0]
            
            new_nodes_1 = len(set(list(df1["GO_term"]))-set(list(df["GO_term"])))
            new_nodes_BP_1 = len(set(list(df1[df1["Ontology_subtype" ] == "biological_process"]["GO_term"]))
                                 -set(list(df[df["Ontology_subtype" ] == "biological_process"]["GO_term"])))
            new_nodes_MF_1 = len(set(list(df1[df1["Ontology_subtype" ] == "molecular_function"]["GO_term"]))
                                 -set(list(df[df["Ontology_subtype" ] == "molecular_function"]["GO_term"])))
            new_nodes_CC_1 = len(set(list(df1[df1["Ontology_subtype" ] == "cellular_component"]["GO_term"]))
                                 -set(list(df[df["Ontology_subtype" ] == "cellular_component"]["GO_term"])))
            old_nodes_1 = len(set(list(df["GO_term"]))-set(list(df1["GO_term"])))
            old_nodes_BP_1 = len(set(list(df[df["Ontology_subtype" ] == "biological_process"]["GO_term"]))
                                 -set(list(df1[df1["Ontology_subtype" ] == "biological_process"]["GO_term"])))
            old_nodes_MF_1 = len(set(list(df[df["Ontology_subtype" ] == "molecular_function"]["GO_term"]))
                                 -set(list(df1[df1["Ontology_subtype" ] == "molecular_function"]["GO_term"])))
            old_nodes_CC_1 = len(set(list(df[df["Ontology_subtype" ] == "cellular_component"]["GO_term"]))
                                 -set(list(df1[df1["Ontology_subtype" ] == "cellular_component"]["GO_term"])))
            obs_df1 = pd.read_csv(path+folder+"_obsolete_nodes",sep = "\t", header= None, names =["GO_term","Ontology_subtype"])
            obs_nodes_1 = obs_df1.shape[0]
            obs_nodes_BP_1 = obs_df1[obs_df1["Ontology_subtype"]=="biological_process"].shape[0]
            obs_nodes_MF_1 = obs_df1[obs_df1["Ontology_subtype"]=="molecular_function"].shape[0]
            obs_nodes_CC_1 = obs_df1[obs_df1["Ontology_subtype"]=="cellular_component"].shape[0]
            new_obs_nodes_1 = len(set(list(obs_df1["GO_term"]))-set(list(obs_df["GO_term"])))
            new_obs_nodes_BP_1 = len(set(list(obs_df1[obs_df1["Ontology_subtype" ] == "biological_process"]["GO_term"]))
                                 -set(list(obs_df[obs_df["Ontology_subtype" ] == "biological_process"]["GO_term"])))
            new_obs_nodes_MF_1 = len(set(list(obs_df1[obs_df1["Ontology_subtype" ] == "molecular_function"]["GO_term"]))
                                 -set(list(obs_df[obs_df["Ontology_subtype" ] == "molecular_function"]["GO_term"])))
            new_obs_nodes_CC_1 = len(set(list(obs_df1[obs_df1["Ontology_subtype" ] == "cellular_component"]["GO_term"]))
                                 -set(list(obs_df[obs_df["Ontology_subtype" ] == "cellular_component"]["GO_term"])))
            old_obs_nodes_1 = len(set(list(obs_df["GO_term"]))-set(list(obs_df1["GO_term"])))
            old_obs_nodes_BP_1 = len(set(list(obs_df[obs_df["Ontology_subtype" ] == "biological_process"]["GO_term"]))
                                 -set(list(obs_df1[obs_df1["Ontology_subtype" ] == "biological_process"]["GO_term"])))

            old_obs_nodes_MF_1 = len(set(list(obs_df[obs_df["Ontology_subtype" ] == "molecular_function"]["GO_term"]))
                                 -set(list(obs_df1[obs_df1["Ontology_subtype" ] == "molecular_function"]["GO_term"])))

            old_obs_nodes_CC_1 = len(set(list(obs_df[obs_df["Ontology_subtype" ] == "cellular_component"]["GO_term"]))
                                 -set(list(obs_df1[obs_df1["Ontology_subtype" ] == "cellular_component"]["GO_term"])))

            

            old_obs_nodes_list = set(list(obs_df["GO_term"]))-set(list(obs_df1["GO_term"]))
            try:
                old_obs_nodes_file = open(path+folder+"_old_obsolete_nodes","w")
            except:
                print("Could not generate old obsoletes file for year" + folder)
                sys.exit()
            for node in old_obs_nodes_list:
                old_obs_nodes_file.write(node+"\n")
            old_obs_nodes_file.close()
            inter = list(set(list(df1["GO_term"])) & set(list(obs_df1["GO_term"])))
            if len(inter) != 0:
                print(folder, "Error")
                sys.exit()
            new_nodes_list = list(set(list(df1["GO_term"]))-set(list(df["GO_term"])))
            for term in new_nodes_list:
                row = {"GO_term":term, "start": int(folder),"end":2024 }
                node_kind = df1[df1["GO_term"] == term]["Ontology_subtype"]
                
                
                if node_kind.shape[0] == 0:
                    continue
                elif node_kind.shape[0] > 1:
                    print(node_kind)
                    print("Error")
                    print(sys.exit())
                else:
                    row["Ontology_subtype"] = node_kind.iloc[0]
                    historical_series.loc[len(historical_series)] = row
                    #print(node_kind.iloc[0])
                    #print(row)

            
            old_nodes_list = list(set(list(df["GO_term"]))-set(list(df1["GO_term"])))
            for term in old_nodes_list:                
                row = historical_series[historical_series["GO_term"] == term]
                if row.shape[0] == 0:
                    continue
                elif row.shape[0] > 1:
                    print(row)
                    print("Error")
                    print(sys.exit())
                else:
                    index = row.index
                    historical_series.iloc[index,3] = int(folder)
                    print(row)

            new_obs_list = set(list(obs_df1["GO_term"]))-set(list(obs_df["GO_term"]))
            for term in new_obs_list:
                row = historical_series[historical_series["GO_term"] == term]
                if row.shape[0] == 0:
                    continue
                elif row.shape[0] > 1:
                    print(row)
                    print("Error")
                    print(sys.exit())
                else:
                    index = row.index
                    historical_series.iloc[index,3] = int(folder)
                    print(row)
            
            print(historical_series)
            data.append([nodes1,nodes_BP_1,nodes_MF_1,nodes_CC_1,
                        new_nodes_1,new_nodes_BP_1,new_nodes_MF_1,new_nodes_CC_1,
                        old_nodes_1,old_nodes_BP_1,old_nodes_MF_1,old_nodes_CC_1,
                        obs_nodes_1,obs_nodes_BP_1,obs_nodes_MF_1,obs_nodes_CC_1,
                        new_obs_nodes_1,new_obs_nodes_BP_1,new_obs_nodes_MF_1,new_obs_nodes_CC_1,
                        old_obs_nodes_1,old_obs_nodes_BP_1,old_obs_nodes_MF_1,old_obs_nodes_CC_1])
            print(data)
            df  = df1
            obs_df = obs_df1


data = pd.DataFrame(data,columns = ["terms","terms_BP","terms_MF","terms_CC",
                                    "new_terms","new_terms_BP","new_terms_MF","new_terms_CC",
                                    "old_terms","old_terms_BP","old_terms_MF","old_terms_CC",
                                    "obsolete_nodes","obsolete_nodes_BP","obsolete_nodes_MF","obsolete_nodes_CC",
                                    "new_obsolete_nodes","new_obsolete_nodes_BP","new_obsolete_nodes_MF","new_obsolete_nodes_CC",
                                    "old_obsolete_nodes","old_obsolete_nodes_BP","old_obsolete_nodes_MF","old_obsolete_nodes_CC"],
                   index= years)

historical_series["term_lifespan"] = historical_series["end"]-historical_series["start"]
#%%
def Z_score(serie):
    return((serie-serie.mean())/serie.std())

def vline_plotter(serie):
    for value in list(serie):
        plt.axvline(value,color="grey", linewidth=1)
def vline_plotter_advance(serie,coord1):
    for value in list(serie):
        axs[coord1].axvline(value,color="grey", linewidth=1.5)
#%%

lifespan_counts = pd.DataFrame(historical_series["term_lifespan"].value_counts())
lifespan_counts["lifespan"] = lifespan_counts.index

lifespan_counts_BP = pd.DataFrame(historical_series[historical_series["Ontology_subtype"]=="biological_process"]["term_lifespan"].value_counts())
lifespan_counts_BP["lifespan"] = lifespan_counts.index
lifespan_counts_BP.sort_values(by="lifespan", ascending = True, inplace = True)

lifespan_counts_MF = pd.DataFrame(historical_series[historical_series["Ontology_subtype"]=="molecular_function"]["term_lifespan"].value_counts())
lifespan_counts_MF["lifespan"] = lifespan_counts.index
lifespan_counts_MF.sort_values(by="lifespan", ascending = True, inplace = True)

lifespan_counts_CC = pd.DataFrame(historical_series[historical_series["Ontology_subtype"]=="cellular_component"]["term_lifespan"].value_counts())
lifespan_counts_CC["lifespan"] = lifespan_counts.index
lifespan_counts_CC.sort_values(by="lifespan", ascending = True, inplace = True)

#Plot historical series
ax = plt.figure(figsize=(12,5))
plt.rc("font",size = 15)
plt.title("Lifespan GO terms")
plt.xlabel("Lifespan (years)")
plt.ylabel("Counts (normalised)")
plt.xticks(range(0,22,1))
ax=plt.plot(lifespan_counts_BP["lifespan"], Z_score(lifespan_counts_BP["count"]), marker = "o", label ="GO:BP",linewidth =  2.5 , color = "blue")
ax=plt.plot(lifespan_counts_MF["lifespan"], Z_score(lifespan_counts_MF["count"]), marker = "o", label ="GO:MF",linewidth =  2.5 , color = "black")
ax=plt.plot(lifespan_counts_CC["lifespan"], Z_score(lifespan_counts_CC["count"]), marker = "o", label ="GO:CC",linewidth =  2.5 , color = "green")
vline_plotter(np.arange(0,22,1))
plt.legend()
plt.tight_layout()
plt.savefig("results/GO_terms_lifespan_Z_score.jpg")

ax = plt.figure(figsize=(12,5))
plt.title("Lifespan GO terms")
plt.rc("font",size = 15)
plt.xlabel("Lifespan (years)")
plt.ylabel("Counts (normalised)")
plt.xticks(range(0,22,1))
ax=plt.plot(lifespan_counts_BP["lifespan"], lifespan_counts_BP["count"]/lifespan_counts_BP["count"].sum(), marker = "o", label ="GO:BP",linewidth =  2.5 , color = "blue")
ax=plt.plot(lifespan_counts_MF["lifespan"], lifespan_counts_MF["count"]/lifespan_counts_MF["count"].sum(), marker = "o", label ="GO:MF",linewidth =  2.5 , color = "black")
ax=plt.plot(lifespan_counts_CC["lifespan"], lifespan_counts_CC["count"]/lifespan_counts_CC["count"].sum(), marker = "o", label ="GO:CC",linewidth =  2.5 , color = "green")
vline_plotter(np.arange(0,22,1))
plt.legend()
plt.tight_layout()
plt.savefig("results/GO_terms_lifespan_normalised.jpg")


#%%

#Plot evolution of terms 
ax = plt.figure(figsize = (12,5))
plt.rc("font",size = 16)
plt.title("Number of GO terms (normalised)")
plt.xlabel("Year")
plt.ylabel("Number of GO terms \n(normalised)")
ax = plt.plot(data.index,Z_score(data["terms_BP"]), marker="o", label ="GO:BP", linewidth =  2.5 , color = "blue")
ax = plt.plot(data.index,Z_score(data["terms_MF"]), marker="o", label ="GO:MF", linewidth =  2.5 ,color = "black")
ax = plt.plot(data.index,Z_score(data["terms_CC"]), marker="o", label ="GO:CC",linewidth =  2.5 , color = "green")
vline_plotter(data.index)
plt.tick_params("x",rotation = 315)
plt.legend()
plt.tight_layout()
#plt.savefig("results/GO_evolution_number_of_terms.svg")
plt.savefig("new_figures_monica/GO_evolution_number_of_terms.jpg", dpi = 300)

"""
#Plot evolution of terms (old)

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of number of GO terms",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of GO terms")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["terms_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["terms_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["terms_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_ylabel("number of terms")


axs[1].title.set_text("Number of GO terms (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["terms_BP"]-data["terms_BP"].mean())/data["terms_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["terms_MF"]-data["terms_MF"].mean())/data["terms_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["terms_CC"]-data["terms_CC"].mean())/data["terms_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_ylabel("number of terms")

plt.savefig("results/GO_evolution_number_of_terms.svg")
"""
#%%

#Plot evolution of new terms and old terms 


ax = plt.figure(figsize = (12,5))
plt.rc("font",size = 16)
#plt.title("Evolution of number of new GO terms (normalised)")
ax = plt.plot(data.index, Z_score(data["new_terms_BP"]),label="GO:BP", marker = "o", linewidth =  2.5 , color = "blue")
ax = plt.plot(data.index, Z_score(data["new_terms_MF"]),label="GO:MF", marker = "o", linewidth =  2.5, color="black")
ax = plt.plot(data.index, Z_score(data["new_terms_CC"]),label="GO:CC", marker = "o",  linewidth =  2.5, color="green")
plt.xlabel("Year")
plt.ylabel("Number of new \n GO terms (normalised)")
vline_plotter(years)
plt.legend()
plt.tick_params("x",rotation = 315)
plt.tight_layout()
#plt.savefig("results/GO_evolution_new_terms.jpg")
plt.savefig("new_figures_monica/GO_evolution_new_terms.svg")

ax = plt.figure(figsize = (12,5))
plt.rc("font",size = 16)
#plt.title("Evolution of number of new obsolete GO terms (normalised)")
ax = plt.plot(data.index, Z_score(data["new_obsolete_nodes_BP"]),label="GO:BP", marker = "o", linewidth =  2.5 , color = "blue")
ax = plt.plot(data.index, Z_score(data["new_obsolete_nodes_MF"]),label="GO:MF", marker = "o", linewidth =  2.5, color="black")
ax = plt.plot(data.index, Z_score(data["new_obsolete_nodes_CC"]),label="GO:CC", marker = "o",  linewidth =  2.5, color="green")
plt.xlabel("Year")
plt.ylabel("Number of new obsolete \n GO terms (normalised)")
vline_plotter(years)
plt.legend()
plt.tick_params("x",rotation = 315)
plt.tight_layout()
#plt.savefig("results/GO_evolution_new_obsolete_terms.jpg")
plt.savefig("new_figures_monica/GO_evolution_new_obsolete_terms.svg")
"""
#Plot evolution of new terms (old)

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(15,15))
plt.suptitle("Changes of GO terms",size = 100)
plt.xlabel("Year")
axs[0].title.set_text("Evolution of number of new GO terms")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_terms_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_terms_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_terms_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_ylabel("number of new terms")
axs[1].title.set_text("Number of new GO terms (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["new_terms_BP"]-data["new_terms_BP"].mean())/data["new_terms_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["new_terms_MF"]-data["new_terms_MF"].mean())/data["new_terms_MF"].std() ,label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["new_terms_CC"]-data["new_terms_CC"].mean())/data["new_terms_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_ylabel("number of new terms")

plt.savefig("results/GO_evolution_number_of_new_terms.svg")



#Plot evolution of old terms

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of number of old GO terms",size = 100)
plt.xlabel("Year")
axs[0].title.set_text("Number of old GO terms")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_terms_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_terms_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_terms_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_ylabel("number of old terms")

axs[1].title.set_text("Number of old GO terms (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["old_terms_BP"]-data["old_terms_BP"].mean())/data["old_terms_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["old_terms_MF"]-data["old_terms_MF"].mean())/data["old_terms_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["old_terms_CC"]-data["old_terms_CC"].mean())/data["old_terms_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_ylabel("number of old terms")


plt.savefig("results/GO_evolution_number_of_old_terms.svg")
"""
#%%

#Plot evolution of obsolete terms

fig, axs  = plt.subplots(3,2, sharex=False, figsize=(60,60))
plt.suptitle("Evolution of obsolete GO terms",size = 100)
plt.xlabel("Year")
axs[0,0].title.set_text("Number of obsolete GO terms")
axs[0,0]=sns.lineplot(data=data, x =data.index, y=data["obsolete_nodes_BP"],label="GO_BP", ax=axs[0,0], linewidth =  2.5 , color = "blue")
axs[0,0]=sns.lineplot(data=data, x =data.index, y=data["obsolete_nodes_MF"],label="GO_MF", ax=axs[0,0], linewidth =  2.5, color="black")
axs[0,0]=sns.lineplot(data=data, x =data.index, y=data["obsolete_nodes_CC"],label="GO_CC", ax=axs[0,0], linewidth =  2.5, color="green")
axs[0,0].set_xlabel("Year")
axs[0,0].set_ylabel("number of obsolete terms")
axs[0,1].title.set_text("Number of obsolete GO terms (normalised)")
axs[0,1]=sns.lineplot(x =data.index, y=(data["obsolete_nodes_BP"]-data["obsolete_nodes_BP"].mean())/data["obsolete_nodes_BP"].std(),label="GO_BP", ax=axs[0,1], linewidth =  2.5 , color = "blue")
axs[0,1]=sns.lineplot(x =data.index, y=(data["obsolete_nodes_MF"]-data["obsolete_nodes_MF"].mean())/data["obsolete_nodes_MF"].std(),label="GO_MF", ax=axs[0,1], linewidth =  2.5, color="black")
axs[0,1]=sns.lineplot(x =data.index, y=(data["obsolete_nodes_CC"]-data["obsolete_nodes_CC"].mean())/data["obsolete_nodes_CC"].std(),label="GO_CC", ax=axs[0,1], linewidth =  2.5, color="green")
axs[0,1].set_xlabel("Year")
axs[0,1].set_ylabel("number of obsolete terms (normalised)")

axs[1,0].title.set_text("Number of new obsolete GO terms")
axs[1,0]=sns.lineplot(data=data, x =data.index, y=data["new_obsolete_nodes_BP"],label="GO_BP", ax=axs[1,0], linewidth =  2.5 , color = "blue")
axs[1,0]=sns.lineplot(data=data, x =data.index, y=data["new_obsolete_nodes_MF"],label="GO_MF", ax=axs[1,0], linewidth =  2.5, color="black")
axs[1,0]=sns.lineplot(data=data, x =data.index, y=data["new_obsolete_nodes_CC"],label="GO_CC", ax=axs[1,0], linewidth =  2.5, color="green")
axs[1,0].set_xlabel("Year")
axs[1,0].set_ylabel("number of new obsolete terms")
axs[1,1].title.set_text("Number of new obsolete GO terms (normalised)")
axs[1,1]=sns.lineplot(x =data.index, y=(data["new_obsolete_nodes_BP"]-data["new_obsolete_nodes_BP"].mean())/data["new_obsolete_nodes_BP"].std(),label="GO_BP", ax=axs[1,1], linewidth =  2.5 , color = "blue")
axs[1,1]=sns.lineplot(x =data.index, y=(data["new_obsolete_nodes_MF"]-data["new_obsolete_nodes_MF"].mean())/data["new_obsolete_nodes_MF"].std(),label="GO_MF", ax=axs[1,1], linewidth =  2.5, color="black")
axs[1,1]=sns.lineplot(x =data.index, y=(data["new_obsolete_nodes_CC"]-data["new_obsolete_nodes_CC"].mean())/data["new_obsolete_nodes_CC"].std(),label="GO_CC", ax=axs[1,1], linewidth =  2.5, color="green")
axs[1,1].set_xlabel("Year")
axs[1,1].set_ylabel("number of new obsolete terms (normalised)")


axs[2,0].title.set_text("Number of old obsolete GO terms")
axs[2,0]=sns.lineplot(data=data, x =data.index, y=data["old_obsolete_nodes_BP"],label="GO_BP", ax=axs[2,0], linewidth =  2.5 , color = "blue")
axs[2,0]=sns.lineplot(data=data, x =data.index, y=data["old_obsolete_nodes_MF"],label="GO_MF", ax=axs[2,0], linewidth =  2.5, color="black")
axs[2,0]=sns.lineplot(data=data, x =data.index, y=data["old_obsolete_nodes_CC"],label="GO_CC", ax=axs[2,0], linewidth =  2.5, color="green")
axs[2,0].set_xlabel("Year")
axs[2,0].set_ylabel("number of old terms")
axs[2,1].title.set_text("Number of old obsolete GO terms (normalised)")
axs[2,1]=sns.lineplot(x =data.index, y=(data["old_obsolete_nodes_BP"]-data["old_obsolete_nodes_BP"].mean())/data["old_obsolete_nodes_BP"].std(),label="GO_BP", ax=axs[2,1], linewidth =  2.5 , color = "blue")
axs[2,1]=sns.lineplot(x =data.index, y=(data["old_obsolete_nodes_MF"]-data["old_obsolete_nodes_MF"].mean())/data["old_obsolete_nodes_MF"].std() ,label="GO_MF", ax=axs[2,1], linewidth =  2.5, color="black")
axs[2,1]=sns.lineplot(x =data.index, y=(data["old_obsolete_nodes_CC"]-data["old_obsolete_nodes_CC"].mean())/data["old_obsolete_nodes_CC"].std(),label="GO_CC", ax=axs[2,1], linewidth =  2.5, color="green")
axs[2,1].set_xlabel("Year")
axs[2,1].set_ylabel("number of old obsolete terms (normalised)")

plt.savefig("results/GO_evolution_obsolete_terms.svg")

#%%

##SAVE HISTORICAL SERIES AND DATA AS CSV

data.to_csv("results/GO_evolution_data_nodes.csv", header = True, index = True)
historical_series.to_csv("results/GO_terms_lifespan.csv", header = True, index = True)












