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
  new_terms = set(list(active_terms_year2)) - set(list(active_terms_year1))
  obsolete_terms = set(list(active_terms_year1)) - set(list(active_terms_year2))
  print(len(new_terms))
  print(len(obsolete_terms))
  print(df[(df["GO_term"].isin(list(obsolete_terms)))& (df["year"]== year1)])


"""
def edge_analyser(graph,edges, year):
    edge_df = pd.DataFrame(edges, columns=['source', 'target', "relationship_type"])
	edge_df
	edge_list = []
    inter_edge_list = []
    for edge in edges:
        edge_kind = edge[2]
        node1 = edge[0]
        node_kind1 = graph.nodes[node1]["namespace"]
        try:
            node_kind1 = conversions[node_kind1]
        except:
            node_kind1 = node_kind1
        node2 = edge[1]
        node_kind2 = graph.nodes[node2]["namespace"]        
        try:
            node_kind2 = conversions[node_kind2]
        except:
            node_kind2 = node_kind2        
        if node_kind1 == node_kind2:
            #print(node1,node2,node_kind1,edge_kind)
            edge_list.append(node1+"\t"+node2+"\t"+node_kind1+"\t"+edge_kind)
        else:
            print("Tuturu. Metaro Upa")
            print(node1,node2,node_kind1,node_kind2,edge_kind)
            inter_edge_list.append(node1+"\t"+node2+"\t"+node_kind1+"\t"+ node_kind2+"\t"+edge_kind)
    return(edge_list,inter_edge_list)
   """     
        

#MAIN RUN

first_year = True
data_folder = os.listdir("/home/input_data/")
print(data_folder)
data_folder.sort()
for file_name in data_folder:
  path = "/home/input_data/"+file_name
  print(path)
  if file_name.endswith (".obo"):
    try:
      graph =  obonet.read_obo(path, ignore_obsolete= False)
    except:
      print("Error opening obo",path)
      sys.exit()
    year = int(file_name.split("_")[0])
    nodes_detailed = node_analyser(graph,year)
#edges_detailed = edge_analyser(graph,graph.edges, year)
    if first_year:
      first_year = False
      all_terms_info = nodes_detailed
      #all_edges_info = edges_detailed
    else:
      all_terms_info= pd.concat([all_terms_info,nodes_detailed])
      node_comparator(all_terms_info, year-1,year)
      break
