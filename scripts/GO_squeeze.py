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

conversions = {"process":"biological_process","function":"molecular_function",
               "component":"cellular_component"}

def node_analyser(graph,nodes, year):
"""    #Three lists: one for nodes, one for obsoletes, one for leaf terms and depth
    node_list = []
    obso_list = []
    leaf_list = []
    first_sons = []"""
    nodes_expanded = []
    for node in nodes:
        node_kind = graph.nodes[node]["namespace"]
        try:
            node_kind = conversions[node_kind]
        except:
            node_kind = node_kind
        if ("is_obsolete" in graph.nodes[node]) and (graph.nodes[node]["is_obsolete"]== "true"):
            is_obs = True
            distance = "None"
        else:
            is_obs = False
            in_deg = graph.in_degree(node)
            if (in_deg == 0) and (graph.degree(node) != 0): #Avoids isolated nodes
                if node_kind == "biological_process":    
                    distance = len(nx.shortest_path(graph, source = node,target = "GO:0008150"))
                elif node_kind == "cellular_component":
                    distance = len(nx.shortest_path(graph, source = node,target = "GO:0005575"))
                elif node_kind == "molecular_function":
                    distance = len(nx.shortest_path(graph, source = node,target = "GO:0003674"))
                else:
                    print("Error", node)
                    sys.exit()
                is_leaf = True
                is_first_layer = False
            elif graph.degree(node) != 0: #first sons
                if node_kind == "biological_process":    
                    distance = len(nx.shortest_path(graph, source = node,target = "GO:0008150"))
                elif node_kind == "cellular_component":
                    distance = len(nx.shortest_path(graph, source = node,target = "GO:0005575"))
                elif node_kind == "molecular_function":
                    distance = len(nx.shortest_path(graph, source = node,target = "GO:0003674"))
                else:
                    print("Error", node)
                    sys.exit()
                  
                if distance == 2:
                  is_leaf = False
                  is_first_layer = True
        nodes_expanded.append([year, node, node_kind, distance, is_leaf, is_first_layer, is_obs])
  #pd.dataframe??
  return(nodes_expanded)

def edge_analyser(graph,edges, year):
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
        
        
def squeezer(file):
    year = file.split("/")[-1].split("_")[0]
    try:
        graph =  obonet.read_obo(file, ignore_obsolete= False)
    except:
        print("Error opening obo",file)
        sys.exit()
    nodes = graph.nodes
    edges = graph.edges
    nodes_detailed,obsoletes,leaf_terms,first_sons = node_analyser(graph,nodes)
    try:
        fh = open(path+"/"+year+"_nodes","w")
    except:
        print("Error opening nodes file:", year)
    fh.write("\n".join(nodes_detailed))
    fh.close()
    try:
        fh = open(path+"/"+year+"_obsolete_nodes","w")
    except:
        print("Error opening nodes file:", year)
    fh.write("\n".join(obsoletes))
    fh.close()
    try:
        fh = open(path+"/"+year+"_leaf_terms","w")
    except:
        print("Error opening nodes file:", year)
    fh.write("\n".join(leaf_terms))
    fh.close()
    
    edges_detailed,inter_edges=edge_analyser(graph,edges)
    try:
        fh = open(path+"/"+year+"_edges","w")
    except:
        print("Error opening edges file:", year)
    fh.write("\n".join(edges_detailed))
    fh.close()
    if len(inter_edges) != 0:
        try:
            fh = open(path+"/"+year+"_inter_edges","w")
        except:
            print("Error opening inter edges file:", year)
        fh.write("\n".join(inter_edges))
        fh.close()
    else:
        print("No inter edges")
    
    try:
        fh = open(path+"/"+year+"_first_sons","w")
    except:
        print("Error opening first sons file:", year)
    fh.write("\n".join(first_sons))
    fh.close()

#MAIN RUN

first_year = True
data_folder = os.listdir("input_data/")
print(data_folder, file=sys.stderr())
data_folder.sort()
for file in data_folder:
    path = "input_data/"+file
    print(path)
    squeezer(path+file)
        year = file.split("/")[-1].split("_")[0]
    try:
        graph =  obonet.read_obo(file, ignore_obsolete= False)
    except:
        print("Error opening obo",file)
        sys.exit()
    nodes_detailed = node_analyser(graph,graph.nodes, year)
    edges_detailed = edge_analyser(graph,graph.edges, year)
    if first_year:
        first_year = False
