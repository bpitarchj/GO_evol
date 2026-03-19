# GO_evol
A pipeline designed to analyse historical evolution in Gene Ontology (GO).


# Libraries required

Complete enviroment used can be checked in .yml file. Specific libraries required are:
 - os
 - sys
 - matplotlib
 - numpy
 - scipy
 - pandas
 - seaborn
 - collections
 - wordcloud
 - itertools
 - obonet
 - networkx

# Input files

The input for analysing the ontology are go.obo and annotation files that can be obtained in the GO Data Archive (https://release.geneontology.org/ ). As all of the files, regardless of the year, are called the same, it is *required to save them adding the year to the file name at the beginning and using "_"*. Example: 2004_go.obo

# Pipeline
## Ontology files (go.obo)
### 1. GO_squeezer.py
Given a folder containing all the go.obo files, GO_squeezer.py will generate a series of text files containing the information of terms, relationships, leaf terms (those at the bottom of the ontology), first layer terms (those directly related to the top of the ontology), obsolete terms and the depth (distance to the top of the ontology) of all the GO terms. The script will one of each of these files per year.
### 2. nodes_plotter.py
Taking the terms and obsolete terms file, this script will generate plots regarding the evolution of the number of terms, number of new terms, new obsolete terms and the lifespan of terms.
### 3. edges_plotter.py
Taking the edges and nodes files, this script will generate plots regarding evolution of relationships and the ratio of relationships/terms over time.
### 4. leaves_plotter.py
Taking the leaf terms information and the depth, this script will generate the ratio of internal_terms/leaf_terms over time, as well as the distribution of depth over time.
### 5. GO_first_sons.py
Taking the first layer terms file, this script will generate the evolution of GO first layer terms over time.
### 6. word_cloud_GO.py
Taking the ontology files directly, this will analyse the new terms and synonyms that appear every year and perform a hypergeometric test to obtain which words are enriched each year.

## Annotation files
### GOA_Uniprot
#### 1. annotation_files_simplifier.py
Only needed for goa_uniprot_files, these script performs a prefiltering step, generating a simplified version of the file with only the columns that are relevant for the analysis.
#### 2. go_annotations_evolution_Uniprot.py
Given the filtered files from previous script, these script will generate plots regarding the evolution in the number of proteins annotated, number of  annotations (with and without IEAs) and the number of annotations by evidence type (with and without IEAs) as described in Gene Ontology (https://geneontology.org/docs/guide-go-evidence-codes/).

### Other sources
#### 1. go_annotations_evolution_MGI.py
Taking the annotation files, this will generate plots regarding the evolution in the number of proteins annotated, number of  annotations and the number of annotations by evidence type as described in Gene Ontology (https://geneontology.org/docs/guide-go-evidence-codes/). The annotation files must be separated in different folders under the name of the source (example: "MGI/" for MGI data) and that source name can be modified to analyse as many organisms as needed.
