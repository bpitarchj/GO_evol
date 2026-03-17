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

org_tag = "MGI"
folder = "/home/bpitarch/Desktop/GO_evolution/annotations_evolutions/annotation_files/"+org_tag+"/"
out_folder= "/home/bpitarch/Desktop/GO_evolution/annotations_evolutions/results/"+org_tag+"/"
#function to loop
files =os.listdir(folder)
files.sort()
first_year = True

for file in files:
    print(file)
    year = file.split("_")[0]
    pre_df = []
    file = folder + file
    with open(file) as f:
        for line in f:
            if not line.startswith("!"):
                line_split = line.rstrip("\n").split("\t")
                #print(line_split)
                gen = line_split[1]
                go_term = line_split[4]
                ann_type = line_split[6]
                go_kind = line_split[8]
                tag = ann_type+"-"+go_kind
                pre_df.append([gen,go_term,ann_type,go_kind,tag])
    df = pd.DataFrame(pre_df)
    proteins = set(df[0])
    annotations_per_go= df[1].value_counts()
    annotations_per_go.name = year
    evidences = df[2].value_counts()
    evidences.name = year
    ontology_subtype= df[3].value_counts()
    ontology_subtype.name = year
    evid_by_onto_type= df[4].value_counts()
    evid_by_onto_type.name = year
    
    df.to_csv(out_folder+year+"_annotations.csv", header= False)
    if first_year:
        global_proteins = [len(proteins)]
        global_annotations_per_go= annotations_per_go
        global_evidences = evidences
        global_ontology_subtype= ontology_subtype
        global_evid_by_onto_type=evid_by_onto_type
        first_year = False
    else:
        global_proteins.append(len(proteins))
        global_annotations_per_go= pd.concat([global_annotations_per_go, annotations_per_go], axis =1)
        global_evidences= pd.concat([global_evidences,evidences],axis =1)
        global_ontology_subtype= pd.concat([global_ontology_subtype,ontology_subtype],axis =1)
        global_evid_by_onto_type= pd.concat([global_evid_by_onto_type,evid_by_onto_type],axis =1)

global_evid_by_onto_type=global_evid_by_onto_type.T
global_evidences = global_evidences.T
global_ontology_subtype = global_ontology_subtype.T
global_annotations_per_go=global_annotations_per_go.T
#sns.lineplot(global_annotations_per_go)# Variation over time of all terms and the distribution
#sns.lineplot(global_ontology_subtype)
#sns.lineplot(global_evid_by_onto_type)
#plt.show()
print(global_proteins)
global_annotations_per_go.to_csv(out_folder+"global_annotations_per_GO.csv")
global_evidences.to_csv(out_folder+"global_evidences.csv")
global_ontology_subtype.to_csv(out_folder+"global_ontology_subtype.csv")
global_evid_by_onto_type.to_csv(out_folder+"global_evid_by_onto_type.csv")
print(global_annotations_per_go)
print(global_evidences)
print(global_ontology_subtype)
print(global_evid_by_onto_type)


#%%
def Z_score(serie):
    return((serie-serie.mean())/serie.std())

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

ax = plt.figure(figsize=(12,5))
plt.title("Number of annotated proteins per year in " + org_tag)
plt.rc("font",size = 15)
plt.xlabel("Year")
plt.ylabel("Number of proteins")
ax = plt.plot(global_evidences.index,global_proteins, marker="o",linewidth=3, color="red" )
plt.tick_params("x",rotation = 315)
vline_plotter(global_evidences.index,"-")
plt.tight_layout()
plt.savefig(out_folder+"annotated_proteins_"+org_tag+".jpg")


"""ax = plt.figure(figsize=(20,12))
plt.title("Number of annotated proteins per year")
plt.xlabel("Year")
plt.ylabel("Number proteins")
ax = plt.bar(global_evidences.index,global_proteins)
plt.savefig(out_folder+"annotated_proteins_"+org_tag+".jpg")"""

#%%

column_groups={"Experimental": ["EXP","IDA","IPI","IMP","IGI","IEP","HDA","HMP"],
               "Phylogeny":["IBA","IKR"],"Computational" : ["ISS","ISO","ISA","ISM","RCA"],
               "Author" : ["TAS","NAS"],"Curator" : ["IC","ND"]}

def column_compactor(df): 
    for group in column_groups.keys():
        columns_to_group = []
        for column in df.columns:
            if column in column_groups[group]:
                columns_to_group.append(column)
        df[group]= df[columns_to_group].sum(axis=1)
        
global_evidences.fillna(0, inplace=True)
column_compactor(global_evidences)
print(global_evidences)
#%%
column_complex_groups={"Experimental-P": ["EXP-P","IDA-P","IPI-P","IMP-P","IGI-P","IEP-P","HDA-P","HMP-P"],
               "Phylogeny-P":["IBA-P","IKR-P"],"Computational-P" : ["ISS-P","ISO-P","ISA-P","ISM-P","RCA-P"],
               "Author-P" : ["TAS-P","NAS-P"],"Curator-P" : ["IC-P","ND-P"],
               "Experimental-F": ["EXP-F","IDA-F","IPI-F","IMP-F","IGI-F","IEP-F","HDA-F","HMP-F"],
               "Phylogeny-F":["IBA-F","IKR-F"],"Computational-F" : ["ISS-F","ISO-F","ISA-F","ISM-F","RCA-F"],
               "Author-F" : ["TAS-F","NAS-F"],"Curator-F" : ["IC-F","ND-F"],
               "Experimental-C": ["EXP-C","IDA-C","IPI-C","IMP-C","IGI-C","IEP-C","HDA-C","HMP-C"],
               "Phylogeny-C":["IBA-C","IKR-C"],"Computational-C" : ["ISS-C","ISO-C","ISA-C","ISM-C","RCA-C"],
               "Author-C" : ["TAS-C","NAS-C"],"Curator-C" : ["IC-C","ND-C"]}
print(column_complex_groups)
def column_complex_compactor(df): 
    for group in column_complex_groups.keys():
        columns_to_group = []
        for column in df.columns:
            if column in column_complex_groups[group]:
                columns_to_group.append(column)
        df[group]= df[columns_to_group].sum(axis=1)

global_evid_by_onto_type.fillna(0, inplace=True)
column_complex_compactor(global_evid_by_onto_type)
print(global_evid_by_onto_type)
#%%
##Number of evidence per ontology subtype

ax = plt.figure(figsize=(12,5))
plt.rc("font",size=16)
#plt.title("Number of annotations in " + org_tag)
plt.title(org_tag)
plt.xlabel("Year")
plt.ylabel("Number of annotations")
ax = plt.plot(global_ontology_subtype.index,global_ontology_subtype["P"], marker="o", label ="GO:BP", linewidth =  2.5 , color = "blue")
ax = plt.plot(global_ontology_subtype.index,global_ontology_subtype["F"], marker="o", label ="GO:MF", linewidth =  2.5 ,color = "black")
ax = plt.plot(global_ontology_subtype.index,global_ontology_subtype["C"], marker="o", label ="GO:CC",linewidth =  2.5 , color = "green")
vline_plotter(global_ontology_subtype.index,"-")
plt.tick_params("x",rotation = 315)
plt.legend()
plt.tight_layout()
#plt.savefig(out_folder+"GO_evolution_annotations_"+org_tag+"_by_ontology_subtype.jpg")
plt.savefig("new_figures_monica/GO_evolution_annotations_"+org_tag+"_by_ontology_subtype.svg")

#%%
#Global evidences Clustering:
    #Experimental = EXP,IDA,IPI,IMP,IGI,IEP,HDA,HMP  (HGI,HEP,HTP are in the GO annotations web, but do not appear in our columns) 
    #Phylo = IBA,IKR (IRD,IBD are in the GO annotations web, but do not appear in our columns) 
    #Computational = ISS,ISO,ISA,ISM,RCA (IGC is in the GO annotations web, but do not appear in our columns) 
    #Author= TAS,NAS
    #Curator = IC, ND
    # IEA
global_evidences.fillna(0, inplace=True)
column_compactor(global_evidences)
print(global_evidences)

ax = plt.figure(figsize=(18,7))
plt.rc("font",size=16)
#plt.title("Number of annotations per evidence type in " + org_tag)
plt.title(org_tag)
plt.xlabel("Year")
plt.ylabel("Number of annotations")
ax = plt.plot(global_evidences.index,global_evidences["Experimental"], marker="o", label ="Experimental", linewidth =  2.5 , color = "blue")
ax = plt.plot(global_evidences.index,global_evidences["Phylogeny"], marker="o", label ="Phylogeny", linewidth =  2.5 ,color = "dimgrey")
ax = plt.plot(global_evidences.index,global_evidences["Computational"], marker="o", label ="Computational",linewidth =  2.5 , color = "turquoise")
ax = plt.plot(global_evidences.index,global_evidences["Author"], marker="o", label ="Author", linewidth =  2.5 ,color = "crimson")
ax = plt.plot(global_evidences.index,global_evidences["Curator"], marker="o", label ="Curator",linewidth =  2.5 , color = "purple")
ax = plt.plot(global_evidences.index,global_evidences["IEA"], marker="o", label ="IEA",linewidth =  2.5 , color = "orange")
vline_plotter(global_evidences.index,"-")
plt.tick_params("x",rotation = 315)
plt.legend()
plt.tight_layout()
#plt.savefig(out_folder+"GO_evolution_clustered_annotations_"+org_tag+".jpg",dpi=250)
plt.savefig("new_figures_monica/GO_evolution_clustered_annotations_"+org_tag+".svg",dpi=250)

#%%
#Global evidences Clustering:
    #Experimental = EXP,IDA,IPI,IMP,IGI,IEP,HDA,HMP  (HGI,HEP,HTP are in the GO annotations web, but do not appear in our columns) 
    #Phylo = IBA,IKR (IRD,IBD are in the GO annotations web, but do not appear in our columns) 
    #Computational = ISS,ISO,ISA,ISM,RCA (IGC is in the GO annotations web, but do not appear in our columns) 
    #Author= TAS,NAS
    #Curator = IC, ND
    # IEA
global_evid_by_onto_type.fillna(0, inplace=True)
column_complex_compactor(global_evid_by_onto_type)
print(global_evid_by_onto_type)

ax = plt.figure(figsize=(16,6))
plt.rc("font",size=15)
plt.title("Number of annotations per evidence type in GO:BP in " + org_tag)
plt.xlabel("Year")
plt.ylabel("Number of annotations")
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Experimental-P"], marker="o", label ="Experimental-P", linewidth =  2.5 , color = "blue")
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Phylogeny-P"], marker="o", label ="Phylogeny-P", linewidth =  2.5 ,color = "dimgrey")
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Computational-P"], marker="o", label ="Computational-P",linewidth =  2.5, color = "turquoise")
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Author-P"], marker="o", label ="Author-P", linewidth =  2.5,color = "crimson" )
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Curator-P"], marker="o", label ="Curator-P",linewidth =  2.5 , color = "purple")
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["IEA-P"], marker="o", label ="IEA-P",linewidth =  2.5 , color = "orange")
vline_plotter(global_evid_by_onto_type.index,"-")
plt.tick_params("x",rotation = 315)
plt.legend()
plt.tight_layout()
plt.savefig(out_folder+"GO_evolution_clustered_annotations_by_GO-BP_"+org_tag+".jpg")

ax = plt.figure(figsize=(16,6))
plt.title("Number of annotations per evidence type in GO:MF in " + org_tag)
plt.xlabel("Year")
plt.ylabel("Number of annotations")
plt.rc("font",size=15)
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Experimental-F"], marker="o", label ="Experimental-F", linewidth =  2.5 , color = "blue")
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Phylogeny-F"], marker="o", label ="Phylogeny-F", linewidth =  2.5,color = "dimgrey" )
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Computational-F"], marker="o", label ="Computational-F",linewidth =  2.5, color = "turquoise" )
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Author-F"], marker="o", label ="Author-F", linewidth =  2.5,color = "crimson" )
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Curator-F"], marker="o", label ="Curator-F",linewidth =  2.5, color = "purple")
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["IEA-F"], marker="o", label ="IEA-F",linewidth =  2.5, color = "orange" )
vline_plotter(global_evid_by_onto_type.index,"-")
plt.tick_params("x",rotation = 315)
plt.legend()
plt.tight_layout()
plt.savefig(out_folder+"GO_evolution_clustered_annotations_by_GO-MF_"+org_tag+".jpg")

ax = plt.figure(figsize=(16,6))
plt.title("Number of annotations per evidence type in GO:CC in " + org_tag)
plt.xlabel("Year")
plt.ylabel("Number of annotations")
plt.rc("font",size=15)
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Experimental-C"], marker="o", label ="Experimental-C", linewidth =  2.5, color = "blue" )
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Phylogeny-C"], marker="o", label ="Phylogeny-C", linewidth =  2.5,color = "dimgrey" )
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Computational-C"], marker="o", label ="Computational-C",linewidth =  2.5, color = "turquoise" )
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Author-C"], marker="o", label ="Author-C", linewidth =  2.5,color = "crimson")
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["Curator-C"], marker="o", label ="Curator-C",linewidth =  2.5, color = "purple" )
ax = plt.plot(global_evid_by_onto_type.index,global_evid_by_onto_type["IEA-C"], marker="o", label ="IEA-C",linewidth =  2.5, color = "orange" )
vline_plotter(global_evid_by_onto_type.index,"-")
plt.tick_params("x",rotation = 315)
plt.legend()
plt.tight_layout()
plt.savefig(out_folder+"GO_evolution_clustered_annotations_by_GO-CC_"+org_tag+".jpg",dpi=250)


#%%

print(global_annotations_per_go)
global_annotations_per_go.fillna(1,inplace = True)
log_global_annotations_per_go = global_annotations_per_go.T.apply(np.log)
print(log_global_annotations_per_go)
ax = plt.figure(figsize=(15,6))
plt.title("Number of annotations per GO term")
plt.xlabel("Years")
plt.ylabel("log10(proteins annotated per GO term)")
sns.boxplot(data= log_global_annotations_per_go)
plt.savefig(out_folder+"GO_evolution_annotations_per_GO_"+org_tag+".jpg")
