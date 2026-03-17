#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:53:45 2024

@author: bpitarch
"""
import os
import sys
import networkx as nx
import obonet
import collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


path = "/home/bpitarch/Desktop/GO_evolution/ontology_evolution/data/"
data_folder = os.listdir(path)
data_folder.sort()
first_year = True

for folder in data_folder:
    print(folder)
    df = pd.read_csv(path+folder+"/"+folder+"_first_sons",sep="\t", header = None)
    df.columns=["GO_term","ontology_subtype"]
    df["year"] = int(folder)
    if first_year:
        fdf = df
        time_scale = df["ontology_subtype"].value_counts()
        first_year=False
    else:
        fdf = pd.concat([fdf,df],axis = 0)
        time_scale = pd.concat([time_scale,df["ontology_subtype"].value_counts()],axis = 1)
print(fdf)
time_scale.columns=data_folder
time_scale = time_scale.T
print(time_scale)
#%%
def Z_score(serie):
    return((serie-serie.mean())/serie.std())

def vline_plotter(serie):
    for value in list(serie):
        plt.axvline(value,color="grey", linewidth=1)
#%%     

#Distribution of direct sons
ax = plt.figure(figsize=(14,5))
#plt.title("Number of terms in the first layer per subontology")
plt.plot(time_scale.index, time_scale["biological_process"],label="GO:BP", marker = "o", linewidth =  2.5 , color = "blue")
plt.plot(time_scale.index,time_scale["molecular_function"] ,label="GO:MF", marker = "o", linewidth =  2.5, color="black")
plt.plot(time_scale.index,time_scale["cellular_component"] ,label="GO:CC", marker = "o", linewidth =  2.5, color="green")
plt.xlabel("Year")
plt.xticks(time_scale.index, rotation=315)
plt.ylabel("number of GO terms")
vline_plotter(time_scale.index)
plt.legend(loc="upper left")
plt.tight_layout()
#plt.savefig("results/GO_first_sons_evolution.svg")
plt.savefig("new_figures_monica/GO_first_sons_evolution.svg")

#%%

terms = list(set(fdf["GO_term"]))

fh = open("results/GO_first_sons_changes","w")
fh.write("GO_term\tOntology_subtype\tYears_as_GO_first_son(comma_separated)\n")

always_first_son_terms = [0,0,0] #BP,MF,CC

for term in terms:
    
    years = list(fdf[fdf["GO_term"]==term]["year"])
    ontology_subtype = fdf[fdf["GO_term"]==term]["ontology_subtype"].iloc[0]
    if len(years) != len(time_scale.index):
        print(term)
        os.system("grep \"" + term + "\" data/*/*obsolete*")
        fh.write(term+"\t"+ontology_subtype+"\t"+",".join([str(year) for year in years])+"\n")
    else:
        print(term+"\t"+ontology_subtype+"\t"+",".join([str(year) for year in years]))
        if ontology_subtype == "biological_process":
            always_first_son_terms[0] +=1
        elif ontology_subtype == "molecular_function":
            always_first_son_terms[1] += 1
        elif ontology_subtype == "cellular_component":
            always_first_son_terms[2] += 1
        else:
            print("Fatal Error")
            sys.exit()
fh.close()
first_sons = np.sum(always_first_son_terms)
print("Number of first sons of all time:",len(terms))
print("Terms that have been always first sons:", first_sons)
print("\tGO:BP constant first sons:", always_first_son_terms[0])
print("\tGO:MF constant first sons:", always_first_son_terms[1])
print("\tGO:CC constant first sons:", always_first_son_terms[2])
##Lifespan??

#%%
"""
descendants = nx.ancestors(graph,"GO:0003674")
print(len(descendants))

for term in descendants:
    
    distance = nx.shortest_path_length(graph,term,"GO:0003674" )
    if (term == "GO:0180024") or (term == "GO:0180020"):
        print(distance)
        print(nx.shortest_path_length(graph,term,term))
    if distance <2:
        if graph.nodes()[term]["namespace"] =="molecular_function":
            #print(term,graph.nodes()[term]["name"])
            depths = []
            for protein in f:
                try:
                    dist = nx.shortest_path_length(graph,protein,term)
                except:
                    continue 
                depths.append(dist)
            if len(depths) > 0:
                print(term,graph.nodes()[term]["name"],max(depths))

"""
               
            
    