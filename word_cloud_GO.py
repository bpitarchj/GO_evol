#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:08:00 2025

@author: bpitarch
"""
import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
import obonet
from collections import Counter
from scipy.stats import hypergeom as hypergeom
import wordcloud
import itertools
conversions = {"process":"biological_process","function":"molecular_function",
               "component":"cellular_component"}

data_folder = os.listdir("data/")
data_folder.sort()
print(data_folder)
years = data_folder
print(years)

subs= [".",",",";",":","(",")","[","]"]

#%%
def word_cleaner(string):
    for pattern in subs:
        string = string.replace(pattern,"")
    return(string.lower())

#%% #GET NEW GO TERMS PER YEAR AND GO CLOUD (normal and by pairs)
# WE NEED PAIRS SEPARATION: HOW DO WE DO IT??? (Term by pairs? What about odd number of words?)
#How about "-" and ":" characters???

new_nodes_per_year = {}

synonym = "exact_synonym"
all_nodes = []
for i in range(0,len(data_folder)):
    folder = data_folder[i]
    path = "data/"+folder+"/"
    file = "data/"+folder+"/"+folder+"_go.obo"
    try:
        graph =  obonet.read_obo(file, ignore_obsolete= True)
    except:
        print("Error opening obo",file)
    nodes_clean = []   
    chiqui_word_cloud = ""
    
    if os.path.isdir(path):
        print(folder)
        with open(path+folder+"_nodes") as f:
            for line in f:
                nodes_clean.append(line.split()[0])
        if i == 0:
            new_nodes_per_year[folder]=nodes_clean
            
            for node in nodes_clean:
                all_nodes.append(node)
                node_info = graph.nodes()[node]["name"]
                chiqui_word_cloud += word_cleaner(node_info) + " "
                try:
                    synonyms = graph.nodes()[node][synonym]
                    for synonym in synonyms:
                        chiqui_word_cloud += word_cleaner(synonym) + " "
                except KeyError:
                    pass
            chiqui_word_cloud = chiqui_word_cloud.split()
            chiqui_pair_word_cloud = [" ".join(chiqui_word_cloud[i:i+2]) for i in range(0, len(chiqui_word_cloud)-1)]
            word_cloud = pd.DataFrame.from_dict(Counter(chiqui_word_cloud), orient="index", columns = [folder])
            pair_word_cloud = pd.DataFrame.from_dict(Counter(chiqui_pair_word_cloud), orient="index", columns = [folder])

        else:
            nodes_cleaner = list(set(nodes_clean)-set(all_nodes))
            new_nodes_per_year[folder]=nodes_cleaner
            for node in nodes_cleaner:  
                all_nodes.append(node)
                node_info = graph.nodes()[node]["name"]
                chiqui_word_cloud += word_cleaner(node_info) + " "
                try:
                    synonyms = graph.nodes()[node][synonym]
                    for synonym in synonyms:
                        chiqui_word_cloud += word_cleaner(synonym) + " "
                except KeyError:
                    pass
            chiqui_word_cloud = chiqui_word_cloud.split()
            chiqui_pair_word_cloud = [" ".join(chiqui_word_cloud[i:i+2]) for i in range(0, len(chiqui_word_cloud)-1)]
            pre_word_cloud = pd.DataFrame.from_dict(Counter(chiqui_word_cloud), orient="index", columns = [folder])
            pre_pair_word_cloud = pd.DataFrame.from_dict(Counter(chiqui_pair_word_cloud), orient="index", columns = [folder])
            word_cloud = pd.concat([word_cloud,pre_word_cloud],axis = 1)
            pair_word_cloud = pd.concat([pair_word_cloud,pre_pair_word_cloud],axis = 1)

            
#%%            
word_cloud.fillna(0.000,inplace = True)
word_cloud["Total"]=word_cloud.sum(axis =1)
word_cloud.loc["Total"]=word_cloud.sum(axis =0)    
pair_word_cloud.fillna(0.000,inplace = True)
pair_word_cloud["Total"]=pair_word_cloud.sum(axis =1)
pair_word_cloud.loc["Total"]=pair_word_cloud.sum(axis =0)      

#%%
word_cloud.to_csv("results/GO_word_cloud.csv",header = True,index= True)
pair_word_cloud.to_csv("results/GO_word_cloud_by pairs.csv",header = True,index= True)

#%%
#In case it colapses, restart from here

word_cloud = pd.read_csv("results/GO_word_cloud.csv",header = 0,index_col= 0)
pair_word_cloud = pd.read_csv("results/GO_word_cloud_by pairs.csv",header = 0,index_col= 0)

print(word_cloud)
print(pair_word_cloud)

#%%
#df['hypergeom_pmf'] = df.apply(lambda row: hypergeom.pmf(row['k'], row['M'], row['n'], row['N']), axis=1)

def hypergeom_calc(df, row, col):
    k = df.loc[row, col]                # word count for that year
    M = df.loc["Total", 'Total']             # total count of all words across all years
    n = df.loc[row, 'Total']          # total count of the word across all years
    N = df.loc['Total', col]            # total count of all words in that year
    return float(hypergeom.sf(k-1, M, n, N))
    #p_val = float(os.popen("Rscript -e '1-phyper("+str(k-1)+","+str(n)+","+str(M-n)+","+ str(N)+")'").read().split()[1])
    #return (p_val)
# Calculate hypergeom for each word-year combination

#%%
result = pd.DataFrame(index=word_cloud.index, columns=word_cloud.columns)
pair_result = pd.DataFrame(index=pair_word_cloud.index, columns=pair_word_cloud.columns)
for row in word_cloud.index:
    for col in word_cloud.columns:
        if row != 'Total' and col != 'Total':
            result.loc[row, col] = hypergeom_calc(word_cloud,row, col)

for row in pair_word_cloud.index:
    for col in pair_word_cloud.columns:
        if row != 'Total' and col != 'Total':
            pair_result.loc[row, col] = hypergeom_calc(pair_word_cloud,row, col)


#%%
result.to_csv("results/hypergeom_stats.csv")
pair_result.to_csv("results/hypergeom_stats_by_pairs.csv")
#%%
#In case it fails
result = pd.read_csv("results/hypergeom_stats.csv", header = 0, index_col = 0)
pair_result = pd.read_csv("results/hypergeom_stats_by_pairs.csv", header = 0, index_col = 0)
#%%
result=result.drop(index="Total", columns="Total")
pair_result=pair_result.drop(index="Total", columns="Total")
#%%

for year in years:
    print(year)
    text = pd.DataFrame(result[result[str(year)]<= 0.001][str(year)], dtype=float)
    text['-log(p_value)'] = np.log10(text+1e-100)*(-1)
    word_freq = dict(zip(text.index,text['-log(p_value)']))

    wc_plot = wordcloud.WordCloud(width=1600, height=800, colormap="inferno", margin=0, mode= "RGBA", background_color=None, max_words=25, relative_scaling=0.2).generate_from_frequencies(word_freq)
    plt.imshow(wc_plot, interpolation='bilinear')
    plt.axis("off")
    plt.title(f'{year} most common words', fontsize=20, color='black', y=1.05)
    plt.margins(x=2, y=8)
    plt.savefig("results/"+str(year)+"_wordcloud.jpg", dpi = 600)
    plt.close()

    text = pd.DataFrame(pair_result[pair_result[str(year)]<= 0.001][str(year)], dtype=float)
    text['-log(p_value)'] = np.log10(text+1e-100)*(-1)
    word_freq = dict(zip(text.index,text['-log(p_value)']))
    wc_plot = wordcloud.WordCloud(width=1600, height=800, colormap="inferno", margin=0, mode= "RGBA", background_color=None, max_words=25, relative_scaling=0.2).generate_from_frequencies(word_freq)
    plt.imshow(wc_plot, interpolation='bilinear')
    plt.title(f'{year} most common pairs of words', fontsize=20, color='black', y=1.05)
    plt.axis("off")
    plt.margins(x=2, y=8)
    plt.savefig("results/"+str(year)+"_wordcloud_by_pairs_of_words.jpg", dpi = 600)
    #print((pair_result[pair_result[year]<= 0.001][year]))
    plt.close()
    
