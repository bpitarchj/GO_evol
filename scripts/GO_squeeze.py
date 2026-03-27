
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 14:32:47 2024

@author: bpitarch
"""
import os
import sys
import networkx as nx
import obonet
import pandas as pd
import matplotlib.pyplot as plt

input_folder = "/home/bpitarch/Desktop/test_go_evol/"

conversions = {"process":"biological_process","function":"molecular_function","component":"cellular_component"}
ontology_top= {"biological_process":"GO:0008150","cellular_component":"GO:0005575","molecular_function":"GO:0003674"}

def node_analyser(graph,year):
  nodes_expanded = []
  for node in graph.nodes:
    #For 2004 GO, names for the subontologies are different. We convert the names to the later code
    node_kind = graph.nodes[node]["namespace"]
    if node_kind in conversions.keys():
      node_kind = conversions[node_kind]
    elif node_kind in conversions.values():
      node_kind = node_kind
    else:
      print("Error: invalid ontology subtype in term:" , node)
    #Check obsoletes
    if ("is_obsolete" in graph.nodes[node]) and (graph.nodes[node]["is_obsolete"]== "true"):
      is_obs = True
      distance = "None"
    else:
      is_obs = False
      if graph.degree(node) == 0:
        print("Isolated term", node)
        #ERROR
      else:
        distance = len(nx.shortest_path(graph, source = node,target = ontology_top[node_kind]))
        if distance == 2:
          is_first_layer = True
        else:
          is_first_layer = False
          in_deg = graph.in_degree(node)
          if (in_deg == 0):
            is_leaf = True
          else:
            is_leaf = False
    nodes_expanded.append([year, node, node_kind, distance, is_leaf, is_first_layer, is_obs])
  nodes_expanded= pd.DataFrame(nodes_expanded, columns=["year","GO_term", "ontology_subtype","distance_to_top","is_leaf","is_first_layer","is_obsolete"])
  #print(nodes_expanded)
  return(nodes_expanded)

def node_comparator(df,year1,year2):
  #Comparison of terms
  active_terms_year1 = df[(df["year"]==year1) & (df["is_obsolete"]==False)]["GO_term"]
  active_terms_year2 = df[(df["year"]==year2) & (df["is_obsolete"]==False)]["GO_term"]
  new_terms = len(set(list(active_terms_year2)) - set(list(active_terms_year1)))
  ####FROM NEW TERMS, GET WORDS FOR WORDCLOUD
  new_obsolete_terms = len(set(list(active_terms_year1)) - set(list(active_terms_year2)))
  return([year2, new_terms,new_obsolete_terms])
    
def edge_squeezer(relationship): #to unify relationship types in "others"
    if relationship == "is_a":
        return(relationship)
    elif relationship == "part_of":
        return(relationship)
    else:
        return("other")
    

def edge_analyser(graph, year): #MAY BE EASILY INTEGRATED IN NODE ANALYSER
  df_columns=["is_a-biological_process","is_a-molecular_function","is_a-cellular_component",
              "part_of-biological_process","part_of-molecular_function","part_of-cellular_component",
              "others-biological_process","others-molecular_function","others-cellular_component",
              "year", "source","interontology", "counts"]
  edge_df = pd.DataFrame(graph.edges, columns=['source', 'target', "relationship_type"])
  node_kinds = nx.get_node_attributes(graph, 'namespace')
  edge_df['source_kind'] = edge_df['source'].map(node_kinds)
  edge_df['target_kind'] = edge_df['target'].map(node_kinds)
  print(edge_df)
  if set(node_kinds.values()) == set(conversions.keys()): #2004
    edge_df['source_kind'] = edge_df['source_kind'].map(conversions)
    edge_df['target_kind'] = edge_df['target_kind'].map(conversions)
  edge_df["interontology"] = edge_df["source_kind"] != edge_df["target_kind"]
  #print(edge_df["interontology"].value_counts())
  #edge_df["Year"] = year
  edge_df["relationship_type_simplified"] = edge_df["relationship_type"].apply(lambda x: edge_squeezer(x))
  
  #final_edge_df = pd.DataFrame(columns=["year","GO_term","ontology_subtype","is_a-BP","is_a-MF","is_a-CC",
  #                                      "part_of-BP","part_of-MF","part_of-CC","others-BP","others-MF", 
  #                                      "others-CC","interontology"])
  edge_df["relationship_ontology_subtype"] = edge_df["relationship_type_simplified"] + "-" + edge_df["source_kind"]
  print(edge_df[edge_df["source_kind"] != edge_df["target_kind"]].groupby("source")["source"].transform("count"))
  sys.exit()
  final_edge_df = edge_df[["source","relationship_ontology_subtype"]].value_counts().reset_index()
  final_edge_df= final_edge_df.pivot(index = "source",columns = "relationship_ontology_subtype", values = "counts")
  final_edge_df = final_edge_df.fillna(0)
  #interontology_counts = edge_df[["source","interontology"]].value_counts().reset_index()
  inter_counts=edge_df[edge_df["source_kind"] != edge_df["target_kind"]].groupby("source")["source"].transform("count")
  print(inter_counts.reindex(edge_df.index).fillna(0).astype(int))
  sys.exit()
  #interontology_counts = edge_df[["source","interontology"]].value_counts().reset_index()
  #interontology_counts = interontology_counts[interontology_counts["interontology"]==True]
  print(interontology_counts)
  #print(final_edge_df.columns)
  #print(final_edge_df)
  final_edge_df = final_edge_df.merge(interontology_counts, left_index=True, right_on="source")
  final_edge_df["year"] = year
  missing_columns = set(df_columns) - set(list(final_edge_df.columns))
  print(missing_columns)
  for column in missing_columns:
    final_edge_df[column] = 0
  return(final_edge_df)



#MAIN RUN

first_year = True
data_folder = os.listdir(input_folder)
print(data_folder)
data_folder.sort()
new_information=[] #Year, New terms and new obsolete terms
for file_name in data_folder:
  path = input_folder+file_name
  print(path)
  if file_name.endswith (".obo"):
    try:
      graph =  obonet.read_obo(path, ignore_obsolete= False)
    except:
      print("Error opening obo",path)
      sys.exit()
    year = int(file_name.split("_")[0])
    nodes_detailed = node_analyser(graph,year)
    edges_detailed = edge_analyser(graph,year)
    continue
    #MERGE NODES AND EDGES DF
    if first_year:
      first_year = False
      all_terms_info = nodes_detailed
      new_terms = len(nodes_detailed[nodes_detailed["is_obsolete"]==False]["GO_term"])
      new_obsolete_terms = len(nodes_detailed[nodes_detailed["is_obsolete"]==True]["GO_term"])
      new_information.append([year,new_terms,new_obsolete_terms])
      #all_edges_info = edges_detailed
    else:
      all_terms_info= pd.concat([all_terms_info,nodes_detailed])
      new_information.append(node_comparator(all_terms_info, year-1,year))

new_terms_info = pd.DataFrame(new_information, columns=["Year","New_terms","New_obsolete_terms"])

#Lifespan

#Plotter

#Plotter functions
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
def leaf_norm(leaves,nodes):
    return((nodes-leaves)/leaves)
