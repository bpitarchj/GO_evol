#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:26:39 2024

@author: bpitarch
"""
#Read as a df. skip first lines
#filter the df, get only protein, term, type of annotation
import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

years = range(2010,2025)
folder = "/home/bpitarch/Desktop/GO_evolution/annotations_evolutions/results/GOA-Uniprot/"

file_types = ["_clustered_evidences.csv","_quality_annotations.csv","_clust_evid_by_GO_type.csv"]

#%%

def vline_plotter(serie, coords):
    if coords == "-": #One normal plot
        for value in list(serie):
            plt.axvline(value,color="grey", linewidth=1)
    elif (type(coords)==list) and (len(coords) == 1):
        for value in list(serie):
            axs[coords[0]].axvline(value,color="grey", linewidth=1.5)
    elif (type(coords)==list) and (len(coords)>1):
            for value in list(serie):
                axs[coords[0],coords[1]].axvline(value,color="grey", linewidth=1.5)
    else:
        print("ERROR")
        sys.exit()
#%%
color_codes={"Experimental": "blue","Phylogeny":"dimgrey","Computational" :"turquoise",
               "Author" : "crimson","Curator" :"purple", "IEA":"orange",
               "Experimental-P": "blue", "Phylogeny-P":"dimgrey","Computational-P" : "turquoise",
               "Author-P" :"crimson" ,"Curator-P" :"purple" ,"IEA-P":"orange",
               "Experimental-F":"blue","Phylogeny-F":"dimgrey","Computational-F" : "turquoise",
               "Author-F" : "crimson","Curator-F" : "purple","IEA-F":"orange",
               "Experimental-C":"blue" ,"Phylogeny-C":"dimgrey","Computational-C" :"turquoise" ,
               "Author-C" : "crimson","Curator-C" : "purple", "IEA-C":"orange"}



#%%
f_type = "_proteins_annotated"
proteins_per_year = [] #First year, then n proteins
swiss_vs_trembl = [] #First year, second n proteins swiss, third n proteins trembl
for year in years:    
    n_proteins = int(os.popen("cat " + folder + str(year)+f_type+"|wc -l").read())
    proteins_per_year.append([year,np.log10(n_proteins)])
    swiss_prots= 0
    trembl_prots = 0
    print(year)
    
    with open(folder + str(year)+f_type,"r") as f:
        for line in f:
            if line.startswith("A0"):
                trembl_prots += 1
            else:
                swiss_prots += 1
    swiss_vs_trembl.append([year,np.log10(swiss_prots),np.log10(trembl_prots)])

protein_array = np.array(proteins_per_year)
swiss_vs_trembl_array = np.array(swiss_vs_trembl)
#%%
ax = plt.figure(figsize=(12,5))
plt.rc("font",size=15)
plt.title("Number of annotated proteins per year")
plt.xlabel("Year")
plt.ylabel("log10 (Number proteins)")
plt.xticks(years,rotation=315)
#plt.bar(protein_array[:,0],protein_array[:,1])
plt.plot(protein_array[:,0],protein_array[:,1], marker="o", color = "red", linewidth =  2.5)
vline_plotter(protein_array[:,0],"-")
plt.tight_layout()
plt.savefig(folder+"annotated_proteins_goa_uniprot.jpg")



ax = plt.figure(figsize=(12,5))
plt.title("Swiss vs TREMBL proteins")
plt.xlabel("Year")
plt.ylabel("log10 (Number proteins)")
plt.xticks(years,rotation =315)
plt.rc("font",size=15)
plt.plot(swiss_vs_trembl_array[:,0],swiss_vs_trembl_array[:,1], marker="o", label ="Swiss-Prot", linewidth =  2.5, color = "gold")
plt.plot(swiss_vs_trembl_array[:,0],swiss_vs_trembl_array[:,2], marker="o", label ="TrEMBL", linewidth =  2.5, color ="grey")
vline_plotter(swiss_vs_trembl_array[:,0],"-")
plt.legend()
plt.tight_layout()
plt.savefig(folder+"swiss_vs_trembl.jpg")



#%%
##Clustered evidences
f_type= "_clustered_evidences.csv"

first_year = True
for year in years:    
    file_name = folder+str(year)+f_type
    df = pd.read_csv(file_name, header = 0, index_col = 0)
    if first_year:
        global_clust_df = df
        first_year = False
    else:
        global_clust_df= pd.concat([global_clust_df, df], axis =1)

global_clust_df.fillna(0, inplace = True)
global_clust_df = global_clust_df.T
global_clust_df.to_csv(folder+"global"+f_type)
print(global_clust_df)

     #all_plots
ax = plt.figure(figsize=(16,6))
plt.rc("font",size = 16)
#plt.title("Number of annotations per evidence type in GOA-Uniprot")
plt.title("GOA-Uniprot")
plt.xlabel("Year")
plt.ylabel("log10(Clustered evidences)")
for column in global_clust_df.columns:
    ax = plt.plot(global_clust_df.index,np.log10(global_clust_df[column]), marker="o", label =column, linewidth =  2.5, color = color_codes[column])
vline_plotter(global_clust_df.index,"-")
plt.legend()
plt.tick_params("x",rotation = 315)
plt.tight_layout()
#plt.savefig(folder+"clustered_evidences_goa_uniprot.png")
plt.savefig("new_figures_monica/clustered_evidences_goa_uniprot.svg")

ax2 = plt.figure(figsize=(12,6))
plt.title("Number of annotations per evidence type in GOA-Uniprot (no IEA)")
plt.xlabel("Year")
plt.ylabel("Clustered evidences")   
for column in global_clust_df.columns:    
 #NO IEA
    if column.startswith("IEA"):
        continue
    else:
        ax2 = plt.plot(global_clust_df.index,global_clust_df[column], marker="o", label =column, linewidth =  2.5, color = color_codes[column])

vline_plotter(global_clust_df.index,"-")
plt.legend()
plt.tight_layout()
plt.tick_params("x",rotation = 315)
plt.savefig(folder+"clustered_evidences_goa_uniprot_no_IEA.png")


#%%

##Clustered evidences by GO
f_type= "_clust_evid_by_GO_type.csv"

first_year = True
for year in years:    
    file_name = folder+str(year)+f_type
    df = pd.read_csv(file_name, header = 0, index_col = 0)
    if first_year:
        global_clust_by_GO_df = df
        first_year = False
    else:
        global_clust_by_GO_df= pd.concat([global_clust_by_GO_df, df], axis =1)

global_clust_by_GO_df.fillna(0, inplace = True)
global_clust_by_GO_df = global_clust_by_GO_df.T
global_clust_by_GO_df.to_csv(folder+"global"+f_type)
print(global_clust_by_GO_df)

     #all_plots
fig1, ax1 = plt.subplots(figsize=(13,6))
ax1.set_title("Number of annotations per evidence type in GO:BP in GOA-Uniprot")
ax1.set_xlabel("Year")
ax1.set_ylabel("log10(Number of annotations)")
fig2, ax2 = plt.subplots(figsize=(13,6))
ax2.set_title("Number of annotations per evidence type in GO:MF in GOA-Uniprot")
ax2.set_xlabel("Year")
ax2.set_ylabel("log10(Number of annotations)")
fig3, ax3 = plt.subplots(figsize=(13,6))
ax3.set_title("Number of annotations per evidence type in GO:CC in GOA-Uniprot")
ax3.set_xlabel("Year")
ax3.set_ylabel("log10(Number of annotations)")

for column in global_clust_by_GO_df.columns:
    if column.endswith("-P"): #Process plot
        ax1.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df[column]), marker="o", label =column, linewidth =  2.5, color = color_codes[column])
    if column.endswith("-F"): #Process plot
        ax2.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df[column]), marker="o", label =column, linewidth =  2.5, color = color_codes[column])
    if column.endswith("-C"): #Process plot
        ax3.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df[column]), marker="o", label =column, linewidth =  2.5, color = color_codes[column])
for year in list(global_clust_by_GO_df.index):    
    ax1.axvline(year,color="grey", linewidth=1)
    ax2.axvline(year,color="grey", linewidth=1)
    ax3.axvline(year,color="grey", linewidth=1)

ax1.tick_params("x",rotation = 315)   
ax2.tick_params("x",rotation = 315) 
ax3.tick_params("x",rotation = 315) 
ax1.legend()
ax2.legend()
ax3.legend()
ax1 = plt.tight_layout()
ax2 = plt.tight_layout()
ax3 = plt.tight_layout()

fig1.savefig(folder+"clustered_evidences_by_GO-BP_goa_uniprot.png")
fig2.savefig(folder+"clustered_evidences_by_GO-MF_goa_uniprot.png")
fig3.savefig(folder+"clustered_evidences_by_GO-CC_goa_uniprot.png")


#NO IEA PLOT
fig4, ax4 = plt.subplots(figsize=(13,6))
ax4.set_title("Number of annotations per evidence type in GO:BP in GOA-Uniprot (no IEA)")
ax4.set_xlabel("Year")
ax4.set_ylabel("log10(Number of annotations)")
fig5, ax5 = plt.subplots(figsize=(13,6))
ax5.set_title("Number of annotations per evidence type in GO:MF in GOA-Uniprot (no IEA)")
ax5.set_xlabel("Year")
ax5.set_ylabel("log10(Number of annotations)")
fig6, ax6 = plt.subplots(figsize=(13,6))
ax6.set_title("Number of annotations per evidence type in GO:CC in GOA-Uniprot (no IEA)")
ax6.set_xlabel("Year")
ax6.set_ylabel("log10(Number of annotations)")

for column in global_clust_by_GO_df.columns:
    #NO IEA
    if column.startswith("IEA"):
           continue
    if column.endswith("-P"): #Process plot
        ax4.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df[column]), marker="o", label =column, linewidth =  2.5, color = color_codes[column])
    elif column.endswith("-F"): #Process plot
        ax5.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df[column]), marker="o", label =column, linewidth =  2.5, color = color_codes[column])
    elif column.endswith("-C"): #Process plot
        ax6.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df[column]), marker="o", label =column, linewidth =  2.5, color = color_codes[column])
for year in list(global_clust_by_GO_df.index):    
    ax4.axvline(year,color="grey", linewidth=1)
    ax5.axvline(year,color="grey", linewidth=1)
    ax6.axvline(year,color="grey", linewidth=1)
    
ax4.tick_params("x",rotation = 315)   
ax5.tick_params("x",rotation = 315) 
ax6.tick_params("x",rotation = 315) 
ax4.legend()
ax5.legend()
ax6.legend()
ax4 = plt.tight_layout()
ax5 = plt.tight_layout()
ax6 = plt.tight_layout()
fig4.savefig(folder+"clustered_evidences_by_GO-BP_goa_uniprot_no_IEA.png")
fig5.savefig(folder+"clustered_evidences_by_GO-MF_goa_uniprot_no_IEA.png")
fig6.savefig(folder+"clustered_evidences_by_GO-CC_goa_uniprot_no_IEA.png")

#%%
global_clust_by_GO_df["GO:BP"] = global_clust_by_GO_df[["Experimental-P","Computational-P","Phylogeny-P","Author-P","Curator-P","IEA-P"]].sum(axis = 1)
global_clust_by_GO_df["GO:MF"] = global_clust_by_GO_df[["Experimental-F","Computational-F","Phylogeny-F","Author-F","Curator-F","IEA-F"]].sum(axis = 1)
global_clust_by_GO_df["GO:CC"] = global_clust_by_GO_df[["Experimental-C","Computational-C","Phylogeny-C","Author-C","Curator-C","IEA-C"]].sum(axis = 1)

fig7, ax7 = plt.subplots(figsize=(12,5))
#ax7.set_title("Evolution of annotations from GOA-Uniprot in GO subontologies")
ax7.set_title("GOA-Uniprot")
plt.rc("font",size = 16)
ax7.set_xlabel("Year")
ax7.set_ylabel("log10(Clustered evidences)")
ax7.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df["GO:BP"]), marker="o", label ="GO:BP", linewidth =  2.5, color = "blue")
ax7.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df["GO:MF"]), marker="o", label ="GO:MF", linewidth =  2.5, color ="black")
ax7.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df["GO:CC"]), marker="o", label ="GO:CC", linewidth =  2.5, color = "green")
for year in list(global_clust_by_GO_df.index):    
    ax7.axvline(year,color="grey", linewidth=1)
ax7.legend()
ax7.tick_params("x",rotation = 315)
ax7 = plt.tight_layout()
#fig7.savefig(folder+"annotations_by_GO_subtype_goa_uniprot.png")
fig7.savefig("new_figures_monica/annotations_by_GO_subtype_goa_uniprot.svg")

global_clust_by_GO_df["GO:BP-no_IEA"] = global_clust_by_GO_df[["Experimental-P","Computational-P","Phylogeny-P","Author-P","Curator-P"]].sum(axis = 1)
global_clust_by_GO_df["GO:MF-no_IEA"] = global_clust_by_GO_df[["Experimental-F","Computational-F","Phylogeny-F","Author-F","Curator-F"]].sum(axis = 1)
global_clust_by_GO_df["GO:CC-no_IEA"] = global_clust_by_GO_df[["Experimental-C","Computational-C","Phylogeny-C","Author-C","Curator-C"]].sum(axis = 1)

fig7, ax7 = plt.subplots(figsize=(12,5))
ax7.set_title("Evolution of annotations from GOA-Uniprot in GO subontologies (no IEA)")
ax7.set_xlabel("Year")
ax7.set_ylabel("log10(Clustered evidences)")
ax7.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df["GO:BP-no_IEA"]), marker="o", label ="GO:BP", linewidth =  2.5, color = "blue")
ax7.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df["GO:MF-no_IEA"]), marker="o", label ="GO:MF", linewidth =  2.5, color ="black")
ax7.plot(global_clust_by_GO_df.index,np.log10(global_clust_by_GO_df["GO:CC-no_IEA"]), marker="o", label ="GO:CC", linewidth =  2.5, color = "green")
for year in list(global_clust_by_GO_df.index):    
    ax7.axvline(year,color="grey", linewidth=1)
ax7.legend()
ax7.tick_params("x",rotation = 315)   
ax7 = plt.tight_layout()
fig7.savefig(folder+"annotations_by_GO_subtype_goa_uniprot_no_IEA.png")

#%%
"""
#Quality annotations
f_type = "_quality_annotations.csv"

first_year = True
for year in years:    
    file_name = folder+str(year)+f_type
    df = pd.read_csv(file_name, header = 0, index_col = 0)
    if first_year:
        global_quality_df = df
        first_year = False
    else:
        global_quality_df= pd.concat([global_quality_df, df], axis =1)

global_quality_df.fillna(0, inplace = True)
global_quality_df = global_quality_df.T
global_quality_df.to_csv(folder+"global"+f_type)
print(global_quality_df)

ax = plt.figure(figsize=(20,12))
plt.title("Source of annotations")
plt.xlabel("Year")
plt.ylabel("Number of annotations")   
for column in global_quality_df.columns:    
    ax = plt.plot(global_quality_df.index,np.log10(global_quality_df[column]), marker="o", label =column, linewidth =  2.5)

vline_plotter(global_quality_df.index,"-")
plt.legend()
plt.savefig(folder+"source_annotations_goa_uniprot.png")
"""
