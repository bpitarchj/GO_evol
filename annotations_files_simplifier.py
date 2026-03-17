import sys
import os
import pandas as pd

input_file = sys.argv[1]
year = input_file.split("/")[-1].split("_")[0]
output_folder=sys.argv[2]
    
pre_df = []
lines = int(os.popen("cat " + input_file + "|wc -l").read().strip())
annotations_per_GO = {}
evidences = {}
ontology_subtype ={}
evid_by_ontology_subtype = {}
uniprot_trembl = {}
clustered_evidences={}
clust_evid_by_GO = {}
proteins = []
counter = 0
#Bigger cluster groups:
#Experimental = EXP,IDA,IPI,IMP,IGI,IEP,HDA,HMP  (HGI,HEP,HTP are in the GO annotations web, but do not appear in our columns) 
#Phylo = IBA,IKR (IRD,IBD are in the GO annotations web, but do not appear in our columns) 
#Computational = ISS,ISO,ISA,ISM,RCA (IGC is in the GO annotations web, but do not appear in our columns) 
#Author= TAS,NAS
#Curator = IC, ND
# IEA

clustered_annotations = {"EXP":"Experimental","IDA":"Experimental","IPI":"Experimental","IMP":"Experimental","IGI":"Experimental","IEP":"Experimental","HDA":"Experimental","HMP":"Experimental","HGI":"Experimental","HEP":"Experimental","HTP":"Experimental","IBA":"Phylogeny","IKR":"Phylogeny","IRD":"Phylogeny","IBD":"Phylogeny","ISS":"Computational","ISO":"Computational","ISA":"Computational","ISM":"Computational","RCA":"Computational","IGC" :"Computational","TAS":"Author","NAS":"Author","IC":"Curator","ND":"Curator","IEA":"IEA"}
def dict_adder(dictionary,key):
	if key in dictionary.keys():
		dictionary[key] += 1
	else:
		dictionary[key] = 1

def dict_closer(dictionary,output_file_name):
	df = pd.DataFrame.from_dict([dictionary]).T
	df.columns = [year]
	df.to_csv(output_folder+"/"+year+"_"+output_file_name)

with open(input_file,"r") as f:
	for line in f:
		counter +=1
		if counter % 1000000 == 0:
			print("Progress:",100*counter/lines, file = sys.stderr)
		line_split = line.rstrip("\n").split("\t")
		quality= line_split[0] #UniprotKB or TrEMBL
		gen = line_split[1]
		go_term = line_split[2]
		ann_type = line_split[3]
		go_kind = line_split[4]
		tag = ann_type+"-"+go_kind
		dict_adder(annotations_per_GO,go_term)
		dict_adder(evidences,ann_type)
		dict_adder(ontology_subtype,go_kind)
		dict_adder(evid_by_ontology_subtype,tag)
		dict_adder(uniprot_trembl,quality)
		clust_ann = clustered_annotations[ann_type]
		dict_adder(clustered_evidences,clust_ann)
		clust_tag = clust_ann+"-"+go_kind
		dict_adder(clust_evid_by_GO,clust_tag)
		proteins.append(gen)

dict_closer(annotations_per_GO, "annotations_per_GO.csv")
dict_closer(evidences, "evidences.csv")
dict_closer(ontology_subtype, "ontology_subtype.csv")
dict_closer(evid_by_ontology_subtype, "evid_by_onto_type.csv")
dict_closer(uniprot_trembl, "quality_annotations.csv")
dict_closer(clustered_evidences, "clustered_evidences.csv")
dict_closer(clust_evid_by_GO, "clust_evid_by_GO_type.csv")
print(year,len(set(proteins)))

