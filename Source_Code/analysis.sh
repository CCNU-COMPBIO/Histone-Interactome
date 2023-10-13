#!/bin/bash
## get uniprot IDs for all histone binding proteins
awk 'BEGIN{FS="\t"}{if($19=="other" && $5!="DNA") print $9}' results_all_histone_PPIs.tsv | sort -u > histone_bp_uniprot_IDs.txt 

##remove duplicates in results
sort -u results_all_histone_PPIs.tsv > results_all.tsv

# Interactome for all organisms

## interactions from nucleosome structures
awk 'BEGIN{FS="\t"}{if($16 == "nucleosome:1") print $0}' results_all.tsv > results_all_nucleosome.tsv
## interactions from histone structures
awk 'BEGIN{FS="\t"}{if(($16 == "nucleosome:0")) print $0}' results_all.tsv > results_all_histone.tsv

## Interactome for human
awk 'BEGIN{FS="\t"}{if($4 == "Homo sapiens") print $0}' results_all.tsv > results_human.tsv
## interactions from nucleosome structures
awk 'BEGIN{FS="\t"}{if($16 == "nucleosome:1") print $0}' results_human.tsv > results_human_nucleosome.tsv
## interactions from histone structures
awk 'BEGIN{FS="\t"}{if(($16 == "nucleosome:0")) print $0}' results_human.tsv >results_human_histone.tsv
