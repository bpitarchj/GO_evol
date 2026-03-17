#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 16:57:43 2024

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
            df = pd.read_csv(path+folder+"_edges",sep = "\t", header= None, names =["term1","term2","Ontology_subtype","relationship_type"])      
            edges = df.shape[0]
            edges_BP = df[df["Ontology_subtype"]=="biological_process"]
            n_edges_BP = edges_BP.shape[0]
            n_edges_BP_isa = edges_BP[edges_BP["relationship_type"]=="is_a"].shape[0]
            n_edges_BP_part = edges_BP[edges_BP["relationship_type"]=="part_of"].shape[0]
            n_edges_BP_others = edges_BP[(edges_BP["relationship_type"]!="is_a") & (edges_BP["relationship_type"]!= "part_of")].shape[0]
            edges_MF = df[df["Ontology_subtype"]=="molecular_function"]
            n_edges_MF = edges_MF.shape[0]
            n_edges_MF_isa = edges_MF[edges_MF["relationship_type"]=="is_a"].shape[0]
            n_edges_MF_part = edges_MF[edges_MF["relationship_type"]=="part_of"].shape[0]
            n_edges_MF_others = edges_MF[(edges_MF["relationship_type"]!="is_a") & (edges_MF["relationship_type"]!= "part_of")].shape[0]
            edges_CC = df[df["Ontology_subtype"]=="cellular_component"]
            n_edges_CC = edges_CC.shape[0]
            n_edges_CC_isa = edges_CC[edges_CC["relationship_type"]=="is_a"].shape[0]
            n_edges_CC_part = edges_CC[edges_CC["relationship_type"]=="part_of"].shape[0]
            n_edges_CC_others = edges_CC[(edges_CC["relationship_type"]!="is_a") & (edges_CC["relationship_type"]!= "part_of")].shape[0]
            n_edges_isa = df[df["relationship_type"]=="is_a"].shape[0]
            n_edges_part = df[df["relationship_type"]=="part_of"].shape[0]
            n_edges_others = df[(df["relationship_type"]!="is_a") & (df["relationship_type"]!= "part_of")].shape[0]

            new_edges = np.nan
            new_edges_isa = np.nan
            new_edges_part = np.nan
            new_edges_others = np.nan
            new_edges_BP = np.nan
            new_edges_BP_isa = np.nan
            new_edges_BP_part = np.nan
            new_edges_BP_others = np.nan
            new_edges_MF = np.nan
            new_edges_MF_isa = np.nan
            new_edges_MF_part = np.nan  
            new_edges_MF_others = np.nan
            new_edges_CC = np.nan
            new_edges_CC_isa = np.nan
            new_edges_CC_part = np.nan
            new_edges_CC_others = np.nan
            
            old_edges = np.nan
            old_edges_isa = np.nan
            old_edges_part = np.nan
            old_edges_others = np.nan
            old_edges_BP = np.nan
            old_edges_BP_isa = np.nan
            old_edges_BP_part = np.nan
            old_edges_BP_others = np.nan           
            old_edges_MF = np.nan
            old_edges_MF_isa = np.nan
            old_edges_MF_part = np.nan
            old_edges_MF_others = np.nan           
            old_edges_CC = np.nan
            old_edges_CC_isa = np.nan
            old_edges_CC_part = np.nan
            old_edges_CC_others = np.nan

            data.append([edges,n_edges_BP,n_edges_MF,n_edges_CC,
                        n_edges_isa,n_edges_BP_isa,n_edges_MF_isa,n_edges_CC_isa,
                        n_edges_part,n_edges_BP_part,n_edges_MF_part,n_edges_CC_part,
                        n_edges_others, n_edges_BP_others,n_edges_MF_others,n_edges_CC_others,
                        new_edges,new_edges_BP,new_edges_MF,new_edges_CC,
                        new_edges_isa,new_edges_BP_isa,new_edges_MF_isa,new_edges_CC_isa,
                        new_edges_part,new_edges_BP_part,new_edges_MF_part,new_edges_CC_part,
                        new_edges_others,new_edges_BP_others,new_edges_MF_others,new_edges_CC_others,
                        old_edges,old_edges_BP,old_edges_MF,old_edges_CC,
                        old_edges_isa,old_edges_BP_isa,old_edges_MF_isa,old_edges_CC_isa,
                        old_edges_part,old_edges_BP_part,old_edges_MF_part,old_edges_CC_part,
                        old_edges_others,old_edges_BP_others,old_edges_MF_others,old_edges_CC_others])
            print(data)
            df["edge"]=df["term1"]+"-"+df["term2"]

        else:
            
            print("Tuturu, metaro upa!")
            df1 = pd.read_csv(path+folder+"_edges",sep = "\t", header= None, names =["term1","term2","Ontology_subtype","relationship_type"])      
            edges1 = df1.shape[0]
            edges_BP1 = df1[df1["Ontology_subtype"]=="biological_process"]
            n_edges_BP1 = edges_BP1.shape[0]
            n_edges_BP_isa1 = edges_BP1[edges_BP1["relationship_type"]=="is_a"].shape[0]
            n_edges_BP_part1 = edges_BP1[edges_BP1["relationship_type"]=="part_of"].shape[0]
            n_edges_BP_others1 = edges_BP1[(edges_BP1["relationship_type"]!="is_a") & (edges_BP1["relationship_type"]!= "part_of")].shape[0]
            edges_MF1 = df1[df1["Ontology_subtype"]=="molecular_function"]
            n_edges_MF1 = edges_MF1.shape[0]
            n_edges_MF_isa1 = edges_MF1[edges_MF1["relationship_type"]=="is_a"].shape[0]
            n_edges_MF_part1 = edges_MF1[edges_MF1["relationship_type"]=="part_of"].shape[0]
            n_edges_MF_others1 = edges_MF1[(edges_MF1["relationship_type"]!="is_a") & (edges_MF1["relationship_type"]!= "part_of")].shape[0]
            edges_CC1 = df1[df1["Ontology_subtype"]=="cellular_component"]
            n_edges_CC1 = edges_CC1.shape[0]
            n_edges_CC_isa1 = edges_CC1[edges_CC1["relationship_type"]=="is_a"].shape[0]
            n_edges_CC_part1 = edges_CC1[edges_CC1["relationship_type"]=="part_of"].shape[0]
            n_edges_CC_others1 = edges_CC1[(edges_CC1["relationship_type"]!="is_a") & (edges_CC1["relationship_type"]!= "part_of")].shape[0]
            n_edges_isa1 = df1[df1["relationship_type"]=="is_a"].shape[0]
            n_edges_part1 = df1[df1["relationship_type"]=="part_of"].shape[0]
            n_edges_others1 = df1[(df1["relationship_type"]!="is_a") & (df1["relationship_type"]!= "part_of")].shape[0]
            
            df1["edge"]=df1["term1"]+"-"+df1["term2"]
            
            new_edges = len(set(list(df1["edge"]))-set(list(df["edge"])))
            new_edges_isa = len(set(list(df1[df1["relationship_type"]=="is_a"]["edge"]))-set(list(df[df["relationship_type"]=="is_a"]["edge"])))
            new_edges_part = len(set(list(df1[df1["relationship_type"]=="part_of"]["edge"]))-set(list(df[df["relationship_type"]=="part_of"]["edge"])))
            new_edges_others =  len(set(list(df1[(df1["relationship_type"]!="is_a") & (df1["relationship_type"]!= "part_of")]["edge"]))-set(list(df[(df["relationship_type"]!="is_a") & (df["relationship_type"]!= "part_of")]["edge"])))
            new_edges_BP = len(set(list(df1[df1["Ontology_subtype"]=="biological_process"]["edge"]))-set(list(df[df["Ontology_subtype"]=="biological_process"]["edge"])))
            new_edges_BP_isa = len(set(list(df1[(df1["Ontology_subtype"]=="biological_process") &(df1["relationship_type"]=="is_a")]["edge"]))-set(list(df[(df["Ontology_subtype"]=="biological_process") & (df["relationship_type"]=="is_a")]["edge"])))
            new_edges_BP_part = len(set(list(df1[(df1["Ontology_subtype"]=="biological_process") &(df1["relationship_type"]=="part_of")]["edge"]))-set(list(df[(df["Ontology_subtype"]=="biological_process") & (df["relationship_type"]=="part_of")]["edge"])))
            new_edges_BP_others = len(set(list(df1[(df1["relationship_type"]!="is_a") & (df1["relationship_type"]!= "part_of") & (df1["Ontology_subtype"]== "biological_process")]["edge"]))-set(list(df[(df["relationship_type"]!="is_a") & (df["relationship_type"]!= "part_of")&(df["Ontology_subtype"]== "biological_process")]["edge"])))
            new_edges_MF = len(set(list(df1[df1["Ontology_subtype"]=="molecular_function"]["edge"]))-set(list(df[df["Ontology_subtype"]=="molecular_function"]["edge"])))
            new_edges_MF_isa = len(set(list(df1[(df1["Ontology_subtype"]=="molecular_function") &(df1["relationship_type"]=="is_a")]["edge"]))-set(list(df[(df["Ontology_subtype"]=="molecular_function") & (df["relationship_type"]=="is_a")]["edge"])))
            new_edges_MF_part = len(set(list(df1[(df1["Ontology_subtype"]=="molecular_function") &(df1["relationship_type"]=="part_of")]["edge"]))-set(list(df[(df["Ontology_subtype"]=="molecular_function") & (df["relationship_type"]=="part_of")]["edge"])))
            new_edges_MF_others = len(set(list(df1[(df1["relationship_type"]!="is_a") & (df1["relationship_type"]!= "part_of") & (df1["Ontology_subtype"]== "molecular_function")]["edge"]))-set(list(df[(df["relationship_type"]!="is_a") & (df["relationship_type"]!= "part_of")&(df["Ontology_subtype"]== "molecular_function")]["edge"])))
            new_edges_CC = len(set(list(df1[df1["Ontology_subtype"]=="cellular_component"]["edge"]))-set(list(df[df["Ontology_subtype"]=="cellular_component"]["edge"])))
            new_edges_CC_isa = len(set(list(df1[(df1["Ontology_subtype"]=="cellular_component") &(df1["relationship_type"]=="is_a")]["edge"]))-set(list(df[(df["Ontology_subtype"]=="cellular_component") & (df["relationship_type"]=="is_a")]["edge"])))
            new_edges_CC_part = len(set(list(df1[(df1["Ontology_subtype"]=="cellular_component") &(df1["relationship_type"]=="part_of")]["edge"]))-set(list(df[(df["Ontology_subtype"]=="cellular_component") & (df["relationship_type"]=="part_of")]["edge"])))
            new_edges_CC_others = len(set(list(df1[(df1["relationship_type"]!="is_a") & (df1["relationship_type"]!= "part_of") & (df1["Ontology_subtype"]== "cellular_component")]["edge"]))-set(list(df[(df["relationship_type"]!="is_a") & (df["relationship_type"]!= "part_of")&(df["Ontology_subtype"]== "cellular_component")]["edge"])))
           
            old_edges = len(set(list(df["edge"]))-set(list(df1["edge"])))
            old_edges_isa = len(set(list(df[df["relationship_type"]=="is_a"]["edge"]))-set(list(df1[df1["relationship_type"]=="is_a"]["edge"])))
            old_edges_part = len(set(list(df[df["relationship_type"]=="part_of"]["edge"]))-set(list(df1[df1["relationship_type"]=="part_of"]["edge"])))
            old_edges_others = len(set(list(df[(df["relationship_type"]!="is_a") & (df["relationship_type"]!= "part_of")]["edge"]))-set(list(df1[(df1["relationship_type"]!="is_a") & (df1["relationship_type"]!= "part_of")]["edge"])))
            old_edges_BP = len(set(list(df[df["Ontology_subtype"]=="biological_process"]["edge"]))-set(list(df1[df1["Ontology_subtype"]=="biological_process"]["edge"])))
            old_edges_BP_isa = len(set(list(df[(df["Ontology_subtype"]=="biological_process") &(df["relationship_type"]=="is_a")]["edge"]))-set(list(df1[(df1["Ontology_subtype"]=="biological_process") & (df1["relationship_type"]=="is_a")]["edge"])))
            old_edges_BP_part = len(set(list(df[(df["Ontology_subtype"]=="biological_process") &(df["relationship_type"]=="part_of")]["edge"]))-set(list(df1[(df1["Ontology_subtype"]=="biological_process") & (df1["relationship_type"]=="part_of")]["edge"])))
            old_edges_BP_others = len(set(list(df[(df["relationship_type"]!="is_a") & (df["relationship_type"]!= "part_of") & (df["Ontology_subtype"]== "biological_process")]["edge"]))-set(list(df1[(df1["relationship_type"]!="is_a") & (df1["relationship_type"]!= "part_of")&(df1["Ontology_subtype"]== "biological_process")]["edge"])))     
            old_edges_MF = len(set(list(df[df["Ontology_subtype"]=="molecular_function"]["edge"]))-set(list(df1[df1["Ontology_subtype"]=="molecular_function"]["edge"])))
            old_edges_MF_isa = len(set(list(df[(df["Ontology_subtype"]=="molecular_function") &(df["relationship_type"]=="is_a")]["edge"]))-set(list(df1[(df1["Ontology_subtype"]=="molecular_function") & (df1["relationship_type"]=="is_a")]["edge"])))
            old_edges_MF_part = len(set(list(df[(df["Ontology_subtype"]=="molecular_function") &(df["relationship_type"]=="part_of")]["edge"]))-set(list(df1[(df1["Ontology_subtype"]=="molecular_function") & (df1["relationship_type"]=="part_of")]["edge"])))
            old_edges_MF_others = len(set(list(df[(df["relationship_type"]!="is_a") & (df["relationship_type"]!= "part_of") & (df["Ontology_subtype"]== "molecular_function")]["edge"]))-set(list(df1[(df1["relationship_type"]!="is_a") & (df1["relationship_type"]!= "part_of")&(df1["Ontology_subtype"]== "molecular_function")]["edge"])))     
            old_edges_CC = len(set(list(df[df["Ontology_subtype"]=="cellular_component"]["edge"]))-set(list(df1[df1["Ontology_subtype"]=="cellular_component"]["edge"])))
            old_edges_CC_isa = len(set(list(df[(df["Ontology_subtype"]=="cellular_component") &(df["relationship_type"]=="is_a")]["edge"]))-set(list(df1[(df1["Ontology_subtype"]=="cellular_component") & (df1["relationship_type"]=="is_a")]["edge"])))
            old_edges_CC_part = len(set(list(df[(df["Ontology_subtype"]=="cellular_component") &(df["relationship_type"]=="part_of")]["edge"]))-set(list(df1[(df1["Ontology_subtype"]=="cellular_component") & (df1["relationship_type"]=="part_of")]["edge"])))
            old_edges_CC_others = len(set(list(df[(df["relationship_type"]!="is_a") & (df["relationship_type"]!= "part_of") & (df["Ontology_subtype"]== "cellular_component")]["edge"]))-set(list(df1[(df1["relationship_type"]!="is_a") & (df1["relationship_type"]!= "part_of")&(df1["Ontology_subtype"]== "cellular_component")]["edge"])))


            data.append([edges1,n_edges_BP1,n_edges_MF1,n_edges_CC1,
                        n_edges_isa1,n_edges_BP_isa1,n_edges_MF_isa1,n_edges_CC_isa1,
                        n_edges_part1,n_edges_BP_part1,n_edges_MF_part1,n_edges_CC_part1,
                        n_edges_others1, n_edges_BP_others1,n_edges_MF_others1,n_edges_CC_others1,
                        new_edges,new_edges_BP,new_edges_MF,new_edges_CC,
                        new_edges_isa,new_edges_BP_isa,new_edges_MF_isa,new_edges_CC_isa,
                        new_edges_part,new_edges_BP_part,new_edges_MF_part,new_edges_CC_part,
                        new_edges_others,new_edges_BP_others,new_edges_MF_others,new_edges_CC_others,
                        old_edges,old_edges_BP,old_edges_MF,old_edges_CC,
                        old_edges_isa,old_edges_BP_isa,old_edges_MF_isa,old_edges_CC_isa,
                        old_edges_part,old_edges_BP_part,old_edges_MF_part,old_edges_CC_part,
                        old_edges_others,old_edges_BP_others,old_edges_MF_others,old_edges_CC_others])
            print(data)
            
            df = df1
        
columns = ["total","BP","MF","CC","is_a","is_a_&_BP","is_a_&_MF","is_a_&_CC",
           "part_of","part_of_&_BP","part_of_&_MF","part_of_&_CC",
           "others","others_&_BP","others_&_MF","others_&_CC",
           "new","new_BP","new_MF","new_CC","new_is_a","new_is_a_&_BP","new_is_a_&_MF","new_is_a_&_CC",
           "new_part_of","new_part_of_&_BP","new_part_of_&_MF","new_part_of_&_CC",
           "new_others","new_others_&_BP","new_others_&_MF","new_others_&_CC",
           "old","old_BP","old_MF","old_CC","old_is_a","old_is_a_&_BP","old_is_a_&_MF","old_is_a_&_CC",
           "old_part_of","old_part_of_&_BP","old_part_of_&_MF","old_part_of_&_CC",
           "old_others","old_others_&_BP","old_others_&_MF","old_others_&_CC",
           ]
data = pd.DataFrame(data, index= years, columns = columns)
print(data)

#%%
def Z_score(serie):
    return((serie-serie.mean())/serie.std())

def vline_plotter(serie):
    for value in list(serie):
        plt.axvline(value,color="grey", linewidth=1)
        
def vline_plotter_medium(serie,coord1):
    for value in list(serie):
        axs[coord1].axvline(value,color="grey", linewidth=1.5)
        
def vline_plotter_advance(serie,coord1,coord2):
    for value in list(serie):
        axs[coord1,coord2].axvline(value,color="grey", linewidth=1.5)
        
def nodes_norm(edges,nodes):
    return(edges/nodes)


nodes = pd.read_csv("results/GO_evolution_data_nodes.csv",header = 0, index_col= 0)
print(nodes)
print(data)
print(data["BP"]/nodes["terms_BP"])
#Normalise everything by number of nodes

#%%
ax = plt.figure(figsize=(12,5))
plt.title("\"is a\" relationships")
plt.rc("font",size = 16)
ax= plt.plot(data.index, Z_score(data["is_a_&_BP"]),label="GO:BP", marker = "o", linewidth =  2.5 , color = "blue")
ax = plt.plot(data.index, Z_score(data["is_a_&_MF"]),label="GO:MF", marker = "o",  linewidth =  2.5, color="black")
ax = plt.plot(data.index, Z_score(data["is_a_&_CC"]),label="GO:CC", marker = "o",  linewidth =  2.5, color="green")
plt.xlabel("Year")
plt.xticks(years, rotation=315)
plt.ylabel("number of relationships \n (normalized)")
plt.tight_layout()
vline_plotter(years)
plt.legend()
#plt.savefig("results/GO_evolution_is_a.jpg")
plt.savefig("new_figures_monica/GO_evolution_is_a.svg")

ax = plt.figure(figsize=(12,5))
plt.rc("font",size = 16)
plt.title("\"is a\" relationships")
ax= plt.plot(data.index,Z_score(nodes_norm(data["is_a_&_BP"], nodes["terms_BP"])),label="GO:BP", marker = "o", linewidth =  2.5 , color = "blue")
ax = plt.plot(data.index, Z_score(nodes_norm(data["is_a_&_MF"], nodes["terms_MF"])),label="GO:MF", marker = "o",  linewidth =  2.5, color="black")
ax = plt.plot(data.index, Z_score(nodes_norm(data["is_a_&_CC"], nodes["terms_CC"])),label="GO:CC", marker = "o",  linewidth =  2.5, color="green")
plt.xlabel("Year")
plt.xticks(years, rotation=315)
plt.ylabel("number of relationships/ \n number of terms (normalized)")
plt.tight_layout()
vline_plotter(years)
plt.legend(loc="upper left")
#plt.savefig("results/GO_evolution_is_a_normalized_terms.jpg")
plt.savefig("new_figures_monica/GO_evolution_is_a_normalized_terms.svg")

ax = plt.figure(figsize=(12,5))
plt.title("\"part of\" relationships")
plt.rc("font",size = 16)
ax= plt.plot(data.index, Z_score(data["part_of_&_BP"]),label="GO:BP", marker = "o", linewidth =  2.5 , color = "blue")
ax = plt.plot(data.index, Z_score(data["part_of_&_MF"]),label="GO:MF", marker = "o",  linewidth =  2.5, color="black")
ax = plt.plot(data.index, Z_score(data["part_of_&_CC"]),label="GO:CC", marker = "o",  linewidth =  2.5, color="green")
plt.xlabel("Year")
plt.xticks(years, rotation=315)
plt.ylabel("number of relationships \n (normalized)")
plt.tight_layout()
vline_plotter(years)
plt.legend()
#plt.savefig("results/GO_evolution_part_of.jpg")
plt.savefig("new_figures_monica/GO_evolution_part_of.svg")

ax = plt.figure(figsize=(12,5))
plt.title("\"part of\" relationships")
plt.rc("font",size = 16)
ax= plt.plot(data.index,Z_score(nodes_norm(data["part_of_&_BP"], nodes["terms_BP"])),label="GO:BP", marker = "o", linewidth =  2.5 , color = "blue")
ax = plt.plot(data.index, Z_score(nodes_norm(data["part_of_&_MF"], nodes["terms_MF"])),label="GO:MF", marker = "o",  linewidth =  2.5, color="black")
ax = plt.plot(data.index, Z_score(nodes_norm(data["part_of_&_CC"], nodes["terms_CC"])),label="GO:CC", marker = "o",  linewidth =  2.5, color="green")
plt.xlabel("Year")
plt.xticks(years, rotation=315)
plt.ylabel("number of relationships/ \n number of terms (normalized)")
plt.tight_layout()
vline_plotter(years)
plt.legend()
#plt.savefig("results/GO_evolution_part_of_normalized_terms.jpg")
plt.savefig("new_figures_monica/GO_evolution_part_of_normalized_terms.svg")

ax = plt.figure(figsize=(12,5))
plt.title("\"others\" relationships")
plt.rc("font",size = 16)
ax= plt.plot(data.index[4:], Z_score(data["others_&_BP"][4:]),label="GO:BP", marker = "o", linewidth =  2.5 , color = "blue")
ax = plt.plot(data.index[14:], Z_score(data["others_&_MF"][14:]),label="GO:MF", marker = "o",  linewidth =  2.5, color="black")
ax = plt.plot(data.index[14:], Z_score(data["others_&_CC"][14:]),label="GO:CC", marker = "o",  linewidth =  2.5, color="green")
plt.xlabel("Year")
plt.xticks(years[4:], rotation=315)
plt.ylabel("number of relationships \n (normalized)")
plt.tight_layout()
vline_plotter(years[4:])
plt.legend()
#plt.savefig("results/GO_evolution_others.jpg")
plt.savefig("new_figures_monica/GO_evolution_others.svg")

ax = plt.figure(figsize=(12,5))
plt.title("\"others\" relationships")
plt.rc("font",size = 16)
ax= plt.plot(data.index[4:],Z_score(nodes_norm(data["others_&_BP"][4:], nodes["terms_BP"][4:])),label="GO:BP", marker = "o", linewidth =  2.5 , color = "blue")
ax = plt.plot(data.index[14:], Z_score(nodes_norm(data["others_&_MF"][14:],nodes["terms_MF"][14:])),label="GO:MF", marker = "o",  linewidth =  2.5, color="black")
ax = plt.plot(data.index[14:], Z_score(nodes_norm(data["others_&_CC"][14:], nodes["terms_CC"][14:])),label="GO:CC", marker = "o",  linewidth =  2.5, color="green")
plt.xlabel("Year")
plt.xticks(years[4:], rotation=315)
plt.ylabel("number of relationships/ \n number of terms (normalized)")
plt.tight_layout()
vline_plotter(years[4:])
plt.legend()
#plt.savefig("results/GO_evolution_others_normalized_terms.jpg")
plt.savefig("new_figures_monica/GO_evolution_others_normalized_terms.svg")

#%%
#Plot evolution of edges

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,25))
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of edges (normalised)")
axs[0].plot(data.index, Z_score(data["BP"]),label="GO_BP", linewidth =  2.5 ,  marker = "o",color = "blue")
axs[0].plot(data.index, Z_score(data["MF"]),label="GO_MF", linewidth =  2.5 ,  marker = "o",color = "black")
axs[0].plot(data.index, Z_score(data["CC"]),label="GO_CC", linewidth =  2.5 ,  marker = "o",color = "green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")
axs[0].legend()
vline_plotter_medium(years,0)


axs[1].title.set_text("Number of edges by number of terms (normalised)")
axs[1].plot(data.index, Z_score(nodes_norm(data["BP"], nodes["terms_BP"])),label="GO_BP", linewidth =  2.5 ,  marker = "o",color = "blue")
axs[1].plot(data.index, Z_score(nodes_norm(data["MF"], nodes["terms_MF"])),label="GO_MF", linewidth =  2.5 ,  marker = "o",color = "black")
axs[1].plot(data.index, Z_score(nodes_norm(data["CC"], nodes["terms_CC"])),label="GO_CC", linewidth =  2.5 ,  marker = "o",color = "green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")
axs[1].legend()
vline_plotter_medium(years,1)

plt.savefig("results/GO_evolution_edges.jpg")
#%%
"""
#%%

#Is a terms

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of \"is a\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of \"is a\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["is_a_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["is_a_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["is_a_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of \"is a\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["is_a_&_BP"]-data["is_a_&_BP"].mean())/data["is_a_&_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["is_a_&_MF"]-data["is_a_&_MF"].mean())/data["is_a_&_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["is_a_&_CC"]-data["is_a_&_CC"].mean())/data["is_a_&_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_is_a_edges.svg")

#%%

#Part of terms

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of \"part of\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of \"part of\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["part_of_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["part_of_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["part_of_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of \"part of\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["part_of_&_BP"]-data["part_of_&_BP"].mean())/data["part_of_&_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y= (data["part_of_&_MF"]-data["part_of_&_MF"].mean())/data["part_of_&_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["part_of_&_CC"]-data["part_of_&_CC"].mean())/data["part_of_&_CC"].std() ,label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_part_of_edges.svg")
#%%
#Others relationships

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of others relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of others relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["others_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["others_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["others_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of others relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["others_&_BP"]-data["others_&_BP"].mean())/data["others_&_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["others_&_MF"]-data["others_&_MF"].mean())/data["others_&_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["others_&_CC"]-data["others_&_CC"].mean())/data["others_&_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_others_edges.svg")

#%%

#Is a terms NEW

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of new \"is a\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of new \"is a\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_is_a_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_is_a_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_is_a_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of new \"is a\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["new_is_a_&_BP"]-data["new_is_a_&_BP"].mean())/data["new_is_a_&_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["new_is_a_&_MF"]-data["new_is_a_&_MF"].mean())/data["new_is_a_&_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["new_is_a_&_CC"]-data["new_is_a_&_CC"].mean())/data["new_is_a_&_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_new_is_a_edges.svg")

#%%
#Part of terms NEW 

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of new \"part of\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of new \"part of\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_part_of_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_part_of_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_part_of_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of new \"part of\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["new_part_of_&_BP"]-data["new_part_of_&_BP"].mean())/data["new_part_of_&_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["new_part_of_&_MF"]-data["new_part_of_&_MF"].mean())/data["new_part_of_&_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["new_part_of_&_CC"]-data["new_part_of_&_CC"].mean())/data["new_part_of_&_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_new_part_of_edges.svg")

#%%
#Others relationships NEW

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of new other relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of new other relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_others_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_others_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_others_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of new other relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["new_others_&_BP"]-data["new_others_&_BP"].mean())/data["new_others_&_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y= (data["new_others_&_MF"]-data["new_others_&_MF"].mean())/data["new_others_&_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["new_others_&_CC"]-data["new_others_&_CC"].mean())/data["new_others_&_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_new_others_edges.svg")

#%%

#Is a terms OLD

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of old \"is a\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of old \"is a\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_is_a_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_is_a_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_is_a_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of old \"is a\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["old_is_a_&_BP"]-data["old_is_a_&_BP"].mean())/data["old_is_a_&_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["old_is_a_&_MF"]-data["old_is_a_&_MF"].mean())/data["old_is_a_&_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["old_is_a_&_CC"]-data["old_is_a_&_CC"].mean())/data["old_is_a_&_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_old_is_a_edges.svg")

#%%

#Part of terms OLD

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of old \"part of\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of old \"part of\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_part_of_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_part_of_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_part_of_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of old \"part of\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["old_part_of_&_BP"]-data["old_part_of_&_BP"].mean())/data["old_part_of_&_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["old_part_of_&_MF"]-data["old_part_of_&_MF"].mean())/data["old_part_of_&_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["old_part_of_&_CC"]-data["old_part_of_&_CC"].mean())/data["old_part_of_&_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_old_part_of_edges.svg")

#%%
#Others relationships OLD

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of old others relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of old others relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_others_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_others_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_others_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of of old others relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=(data["old_others_&_BP"]-data["old_others_&_BP"].mean())/data["old_others_&_BP"].std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=(data["old_others_&_MF"]-data["old_others_&_MF"].mean())/data["old_others_&_MF"].std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=(data["old_others_&_CC"]-data["old_others_&_CC"].mean())/data["old_others_&_CC"].std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_old_others_edges.svg")

#%%

nodes = pd.read_csv("results/GO_evolution_data_nodes.csv",header = 0, index_col= 0)
print(nodes)
print(data)
print(data["BP"]/nodes["terms_BP"])
#Normalise everything by number of nodes
#%%

#Plot evolution of edges by nodes

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of number of GO edges",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of edges")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["BP"]/nodes["terms_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["MF"]/nodes["terms_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["CC"]/nodes["terms_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of edges (normalised)")
axs[1]=sns.lineplot(x =data.index, y=((data["BP"]/nodes["terms_BP"])-data["BP"]/nodes["terms_BP"].mean())/(data["BP"]/nodes["terms_BP"]).std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=((data["MF"]/nodes["terms_MF"])-data["MF"]/nodes["terms_MF"].mean())/(data["MF"]/nodes["terms_MF"]).std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=((data["CC"]/nodes["terms_CC"])-data["CC"]/nodes["terms_CC"].mean())/(data["CC"]/nodes["terms_CC"]).std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_edges_by_nodes.svg")

#%%

#Is a terms by nodes

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of \"is a\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of \"is a\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["is_a_&_BP"]/nodes["terms_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["is_a_&_MF"]/nodes["terms_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["is_a_&_CC"]/nodes["terms_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of \"is a\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=((data["is_a_&_BP"]/nodes["terms_BP"])-data["is_a_&_BP"]/nodes["terms_BP"].mean())/(data["is_a_&_BP"]/nodes["terms_BP"]).std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=((data["is_a_&_MF"]/nodes["terms_MF"])-data["is_a_&_MF"]/nodes["terms_MF"].mean())/(data["is_a_&_MF"]/nodes["terms_MF"]).std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=((data["is_a_&_CC"]/nodes["terms_CC"])-data["is_a_&_CC"]/nodes["terms_CC"].mean())/(data["is_a_&_CC"]/nodes["terms_CC"]).std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_is_a_edges_by_nodes.svg")

#%%

#Part of terms by nodes

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of \"part of\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of \"part of\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["part_of_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["part_of_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["part_of_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of \"part of\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=((data["part_of_&_BP"]/nodes["terms_BP"])-data["part_of_&_BP"]/nodes["terms_BP"].mean())/(data["part_of_&_BP"]/nodes["terms_BP"]).std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=((data["part_of_&_MF"]/nodes["terms_MF"])-data["part_of_&_MF"]/nodes["terms_MF"].mean())/(data["part_of_&_MF"]/nodes["terms_MF"]).std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=((data["part_of_&_CC"]/nodes["terms_CC"])-data["part_of_&_CC"]/nodes["terms_CC"].mean())/(data["part_of_&_CC"]/nodes["terms_CC"]).std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_part_of_edges_by_nodes.svg")
#%%
#Others relationships

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of others relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of others relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["others_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["others_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["others_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of others relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=((data["others_&_BP"]/nodes["terms_BP"])-data["others_&_BP"]/nodes["terms_BP"].mean())/(data["others_&_BP"]/nodes["terms_BP"]).std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=((data["others_&_MF"]/nodes["terms_MF"])-data["others_&_MF"]/nodes["terms_MF"].mean())/(data["others_&_MF"]/nodes["terms_MF"]).std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=((data["others_&_CC"]/nodes["terms_CC"])-data["others_&_CC"]/nodes["terms_CC"].mean())/(data["others_&_CC"]/nodes["terms_CC"]).std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_others_edges_by_nodes.svg")

#%%

#Is a terms NEW by nodes

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of new \"is a\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of new \"is a\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_is_a_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_is_a_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_is_a_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of new \"is a\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=((data["new_is_a_&_BP"]/nodes["terms_BP"])-data["new_is_a_&_BP"]/nodes["terms_BP"].mean())/(data["new_is_a_&_BP"]/nodes["terms_BP"]).std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=((data["new_is_a_&_MF"]/nodes["terms_MF"])-data["new_is_a_&_MF"]/nodes["terms_MF"].mean())/(data["new_is_a_&_MF"]/nodes["terms_MF"]).std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=((data["new_is_a_&_CC"]/nodes["terms_CC"])-data["new_is_a_&_CC"]/nodes["terms_CC"].mean())/(data["new_is_a_&_CC"]/nodes["terms_CC"]).std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_new_is_a_edges_by_nodes.svg")

#%%
#Part of terms NEW by nodes

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of new \"part of\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of new \"part of\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_part_of_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_part_of_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_part_of_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of new \"part of\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=((data["new_part_of_&_BP"]/nodes["terms_BP"])-data["new_part_of_&_BP"]/nodes["terms_BP"].mean())/(data["new_part_of_&_BP"]/nodes["terms_BP"]).std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=((data["new_part_of_&_MF"]/nodes["terms_MF"])-data["new_part_of_&_MF"]/nodes["terms_MF"].mean())/(data["new_part_of_&_MF"]/nodes["terms_MF"]).std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=((data["new_part_of_&_CC"]/nodes["terms_CC"])-data["new_part_of_&_CC"]/nodes["terms_CC"].mean())/(data["new_part_of_&_CC"]/nodes["terms_CC"]).std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_new_part_of_edges_by_nodes.svg")

#%%
#Others relationships NEW by nodes

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of new other relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of new other relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_others_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_others_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["new_others_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of new other relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=((data["new_others_&_BP"]/nodes["terms_BP"])-data["new_others_&_BP"]/nodes["terms_BP"].mean())/(data["new_others_&_BP"]/nodes["terms_BP"]).std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=((data["new_others_&_MF"]/nodes["terms_MF"])-data["new_others_&_MF"]/nodes["terms_MF"].mean())/(data["new_others_&_MF"]/nodes["terms_MF"]).std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=((data["new_others_&_CC"]/nodes["terms_CC"])-data["new_others_&_CC"]/nodes["terms_CC"].mean())/(data["new_others_&_CC"]/nodes["terms_CC"]).std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_new_others_edges_by_nodes.svg")

#%%

#Is a terms OLD by nodes

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of old \"is a\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of old \"is a\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_is_a_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_is_a_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_is_a_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of old \"is a\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=((data["old_is_a_&_BP"]/nodes["terms_BP"])-data["old_is_a_&_BP"]/nodes["terms_BP"].mean())/(data["old_is_a_&_BP"]/nodes["terms_BP"]).std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=((data["old_is_a_&_MF"]/nodes["terms_MF"])-data["old_is_a_&_MF"]/nodes["terms_MF"].mean())/(data["old_is_a_&_MF"]/nodes["terms_MF"]).std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=((data["old_is_a_&_CC"]/nodes["terms_CC"])-data["old_is_a_&_CC"]/nodes["terms_CC"].mean())/(data["old_is_a_&_CC"]/nodes["terms_CC"]).std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_old_is_a_edges_by_nodes.svg")

#%%

#Part of terms OLD by nodes

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of old \"part of\" relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of old \"part of\" relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_part_of_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_part_of_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_part_of_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of old \"part of\" relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=((data["old_part_of_&_BP"]/nodes["terms_BP"])-data["old_part_of_&_BP"]/nodes["terms_BP"].mean())/(data["old_part_of_&_BP"]/nodes["terms_BP"]).std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=((data["old_part_of_&_MF"]/nodes["terms_MF"])-data["old_part_of_&_MF"]/nodes["terms_MF"].mean())/(data["old_part_of_&_MF"]/nodes["terms_MF"]).std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=((data["old_part_of_&_CC"]/nodes["terms_CC"])-data["old_part_of_&_CC"]/nodes["terms_CC"].mean())/(data["old_part_of_&_CC"]/nodes["terms_CC"]).std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_old_part_of_edges_by_nodes.svg")

#%%
#Others relationships OLD by nodes

fig, axs  = plt.subplots(2,1, sharex=False, figsize=(40,40))
plt.suptitle("Evolution of old others relationships",size = 100)
axs[0].set_xlabel("Year")

axs[0].title.set_text("Number of old others relationships")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_others_&_BP"],label="GO_BP", ax=axs[0], linewidth =  2.5 , color = "blue")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_others_&_MF"],label="GO_MF", ax=axs[0], linewidth =  2.5, color="black")
axs[0]=sns.lineplot(data=data, x =data.index, y=data["old_others_&_CC"],label="GO_CC", ax=axs[0], linewidth =  2.5, color="green")
axs[0].set_xlabel("Year")
axs[0].set_xticks(years)
axs[0].set_ylabel("number of edges")


axs[1].title.set_text("Number of of old others relationships (normalised)")
axs[1]=sns.lineplot(x =data.index, y=((data["old_others_&_BP"]/nodes["terms_BP"])-data["old_others_&_BP"]/nodes["terms_BP"].mean())/(data["old_others_&_BP"]/nodes["terms_BP"]).std(),label="GO_BP", ax=axs[1], linewidth =  2.5 , color = "blue")
axs[1]=sns.lineplot(x =data.index, y=((data["old_others_&_MF"]/nodes["terms_MF"])-data["old_others_&_MF"]/nodes["terms_MF"].mean())/(data["old_others_&_MF"]/nodes["terms_MF"]).std(),label="GO_MF", ax=axs[1], linewidth =  2.5, color="black")
axs[1]=sns.lineplot(x =data.index, y=((data["old_others_&_CC"]/nodes["terms_CC"])-data["old_others_&_CC"]/nodes["terms_CC"].mean())/(data["old_others_&_CC"]/nodes["terms_CC"]).std(),label="GO_CC", ax=axs[1], linewidth =  2.5, color="green")
axs[1].set_xlabel("Year")
axs[1].set_xticks(years)
axs[1].set_ylabel("number of edges")

plt.savefig("results/GO_evolution_old_others_edges_by_nodes.svg")
"""