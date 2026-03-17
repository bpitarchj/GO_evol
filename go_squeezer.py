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
def node_analyser(graph,nodes):
    #Three lists: one for nodes, one for obsoletes, one for leaf terms and depth
    node_list = []
    obso_list = []
    leaf_list = []
    first_sons = []
    for node in nodes:
        node_kind = graph.nodes[node]["namespace"]
        try:
            node_kind = conversions[node_kind]
        except:
            node_kind = node_kind
        if ("is_obsolete" in graph.nodes[node]) and (graph.nodes[node]["is_obsolete"]== "true"):
            #print(node,node_kind,"Obsolete") #Append to obsoletes
            obso_list.append(node+"\t"+node_kind)
        else:
            #print(node,node_kind) #Append to nodes
            node_list.append(node+"\t"+node_kind)
            deg = graph.in_degree(node)
            if (deg == 0) and (graph.degree(node) != 0): #Avoids isolated nodes
                if node_kind == "biological_process":    
                    distance = len(nx.shortest_path(graph, source = node,target = "GO:0008150"))
                elif node_kind == "cellular_component":
                    distance = len(nx.shortest_path(graph, source = node,target = "GO:0005575"))
                elif node_kind == "molecular_function":
                    distance = len(nx.shortest_path(graph, source = node,target = "GO:0003674"))
                else:
                    print("Error", node)
                    sys.exit()
                #print(node,node_kind,distance)#Append
                leaf_list.append(node+"\t"+node_kind+"\t"+str(distance))
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
                    first_sons.append(node+"\t"+node_kind)
    return(node_list,obso_list,leaf_list,first_sons)

def edge_analyser(graph,edges):
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
    path = "/".join(file.split("/")[:-1])
    year = file.split("/")[1]
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
    
data_folder = os.listdir("data/")
print(data_folder)
data_folder.sort()
for folder in data_folder:
    path = "data/"+folder+"/"
    if os.path.isdir(path):
        file = folder +"_go.obo"
        print(file)
        squeezer(path+file)