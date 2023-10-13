#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 10:35:51 2019

@author: pengy10
"""
import urllib
import requests
from time import sleep
import re
import json
#from rcsbsearch import TextQuery
import time
import pandas as pd


##rcsbsearch package is not working, beacasue RCSB PDB Search API has been updated to V2. 
#def search_text(keyword,pdb_list):
#    results = TextQuery(keyword) \
#    .and_("rcsb_entry_info.polymer_entity_count_protein").greater_or_equal(2) \
#    .exec("entry")

#    for entry in results:
#        pdb_list.append(entry)
        

def search_text(keyword,pdb_list):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    search_request = {
  "query": {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
     {
    "type": "terminal",
    "service": "full_text",
    "parameters": {
      "value": keyword
      }
    },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
        "operator": "greater_or_equal",
        "value": 2
        }
      }
    ]
  },
  "request_options": {
   "return_all_hits": True
  },
  "return_type": "entry"
}
        
    response = requests.post(url, json=search_request)
    
    if response.status_code == 200:
    #    print(json.loads(response.text) )
        output= json.loads(response.text)
        for element in output['result_set']:
    #        if(element != ""):
            pdb_list.append(element['identifier'])
    else:
        print("error:404")
        
        
def get_entityIDs(pdb_id,EntityIDs_list):
    url = "https://data.rcsb.org/graphql?query="
    
    search_request = """
    {
       entries(entry_ids: \"""" +str(pdb_id)+ """\") { 
       rcsb_entry_container_identifiers{
        entity_ids
        }
     }
     }
    """
    
    query  = urllib.parse.quote_plus(search_request)
    
    response = requests.get(url+query)
    if response.status_code == 200:
    #    print(json.loads(response.text) )
        output= json.loads(response.text)
        entity_ids = output["data"]['entries'][0]["rcsb_entry_container_identifiers"]['entity_ids']
        for element in entity_ids:
            EntityIDs_list.append(str(pdb_id) + "_" + str(element))        
    else:
        print("error:404")


def generate_customReport(entity,chains_report):
    url = "https://data.rcsb.org/graphql?query="
#    entity = "6O7G_2"
    search_request = """
    {
      polymer_entities(entity_ids: \"""" +str(entity)+ """\") { 
       rcsb_polymer_entity{
        pdbx_description
        }
        entity_poly {
          pdbx_strand_id
          rcsb_entity_polymer_type
        }
        rcsb_entity_source_organism {
          ncbi_scientific_name
        }
        rcsb_polymer_entity_container_identifiers{
        entry_id
        uniprot_ids
        }
       uniprots{
        rcsb_uniprot_protein{
          name
          {value}
        source_organism	{
        scientific_name
          }
        }
    }
      }
     }
    """
    
    query  = urllib.parse.quote_plus(search_request)
    
    response = requests.get(url+query)
    
    if response.status_code == 200:
    #    print(json.loads(response.text) )
        output= json.loads(response.text)
        
        for element in output["data"]['polymer_entities']:
           chain_ID = element['entity_poly']['pdbx_strand_id'].split(",")
           PDB_ID = element['rcsb_polymer_entity_container_identifiers']['entry_id']
           
           if(element['uniprots']):
               uniprot_name = element['uniprots'][0]['rcsb_uniprot_protein']['name']['value']
               source = element['uniprots'][0]['rcsb_uniprot_protein']['source_organism']['scientific_name']
           elif(element['rcsb_entity_source_organism']): 
               source = element['rcsb_entity_source_organism'][0]['ncbi_scientific_name']
               uniprot_name = "NA"
           else:    
               uniprot_name = "NA"
               source = "NA"
              
               
           if(element['rcsb_polymer_entity_container_identifiers']['uniprot_ids']):
               uniprotAcc = ';'.join(element['rcsb_polymer_entity_container_identifiers']['uniprot_ids'])    
           else:
               uniprotAcc = "NA"   
               
           if(element['rcsb_polymer_entity']['pdbx_description']):
               pdbx_description = element['rcsb_polymer_entity']['pdbx_description']
           else:
               pdbx_description = "NA"  
               
           if(element['entity_poly']['rcsb_entity_polymer_type']):
               entityMacromoleculeType = element['entity_poly']['rcsb_entity_polymer_type']
           else:
               entityMacromoleculeType = "NA"           
           
    #       print(chain_ID)
           for ID in chain_ID:
               chains_report.append(str(PDB_ID) + "\t" + str(ID) + "\t" +\
                   str(pdbx_description) + "\t" + str(source) + "\t" + str(uniprotAcc) + "\t" +\
                       str(uniprot_name) + "\t" + str(entityMacromoleculeType))
           
    else:
        print ("error:404")



def search_uniprot(uniprot_ID,pdb_list):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    search_request = {
  "query": {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": uniprot_ID,
          "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
        }
      },
  "request_options": {
    "return_all_hits": True
  },
  "return_type": "entry"
}
        
    response = requests.post(url, json=search_request)
    
    if response.status_code == 200:
    #    print(json.loads(response.text) )
        output= json.loads(response.text)
        for element in output['result_set']:
    #        if(element != ""):
            pdb_list.append(element['identifier'])
    else:
        print("error:404")


def get_CDD_domain_hit(chain_file):
#    with open(chain_file, "r") as fh:
    df = pd.read_table(chain_file, sep="\t")
    num_per_job = 50 # number of PDB chains per job
    job_num = len(df)//num_per_job # totoal number of jobs to be submitted
    CDD_report = list()
    failed_case = list()
    
    for i in range(0,job_num):
        print(i)
        pdb_chain_list = list()
        for j in range(i*num_per_job,(i+1)*num_per_job):
            PDB = df.iloc[j,0]
            chain = df.iloc[j,1]
            pdb_chain_list.append(str(PDB) + "_"  + str(chain))
        #generate query chains for submiting job 
        query_chains = "%0A".join(pdb_chain_list)
        query_text ="queries="+ query_chains + "&tdata=hits"
        urlCustom = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?"
        response1 = requests.post(urlCustom, data=query_text)
        if(response1.status_code == 200):
            try:
                job_id = response1.text.split('\n')[1].split('\t')[1]
                urlCustom = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?cdsid=" + job_id
                response2 = requests.get(urlCustom)                    
                status_code = int(response2.text.split('\n')[3].split('\t')[1])
                time_start = time.time()
                time_duration = 600
                while status_code ==3: ## job is still running:
                    sleep(1)
                    response2 = requests.get(urlCustom)
                    status_code = int(response2.text.split('\n')[3].split('\t')[1])
                    if time.time() > time_start + time_duration:
                        status_code = 6
                if(status_code ==0): ## Job is done successfully
                    if(response2.text.split('\n')[8:] != [""]): ## domain hits available from output
                        CDD_report.append(response2.text.split('\n')[8:])
    #                    print(CDD_report)
                elif(status_code) ==1: ## Invalid search ID
                    failed_case.append("Invalid search ID")
                elif(status_code) ==2: ## No effective input (usually no query proteins or search ID specified)
                    failed_case.append("No effective input")
                elif(status_code) ==4: ## Queue manager (qman) service error
                    failed_case.append("Queue manager (qman) service error")
                elif(status_code) ==5: ## Data is corrupted or no longer available (cache cleaned, etc)
                    failed_case.append("Data is corrupted or no longer available")
                elif(status_code) ==6: ## Data is corrupted or no longer available (cache cleaned, etc)
                    print("job running for too long time")
                    failed_case.append("job running for too long time")
            except:
                print("Something else went wrong in CDD query")
                continue
    ## loop over the rest jobs    
    pdb_chain_list = list()
    for j in range(job_num*num_per_job,len(df)):
        PDB = df.iloc[j,0]
        chain = df.iloc[j,1]
        pdb_chain_list.append(str(PDB) + "_"  + str(chain))
    query_chains = "%0A".join(pdb_chain_list)
    query_text ="queries="+ query_chains + "&tdata=hits"
    urlCustom = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?"
    response1 = requests.post(urlCustom, data=query_text)
    if(response1.status_code == 200):
        try:
            job_id = response1.text.split('\n')[1].split('\t')[1]
            urlCustom = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?cdsid=" + job_id
            response2 = requests.get(urlCustom)                    
            status_code = int(response2.text.split('\n')[3].split('\t')[1])
            time_start = time.time()
            time_duration = 600
            while status_code ==3: ## job is still running:
                sleep(1)
                response2 = requests.get(urlCustom)
                status_code = int(response2.text.split('\n')[3].split('\t')[1])
                if time.time() > time_start + time_duration:
                    status_code = 6
            if(status_code ==0): ## Job is done successfully
                if(response2.text.split('\n')[8:] != [""]): ## domain hits available from output
                    CDD_report.append(response2.text.split('\n')[8:])
#                    print(CDD_report)
            elif(status_code) ==1: ## Invalid search ID
                failed_case.append("Invalid search ID")
            elif(status_code) ==2: ## No effective input (usually no query proteins or search ID specified)
                failed_case.append("No effective input")
            elif(status_code) ==4: ## Queue manager (qman) service error
                failed_case.append("Queue manager (qman) service error")
            elif(status_code) ==5: ## Data is corrupted or no longer available (cache cleaned, etc)
                failed_case.append("Data is corrupted or no longer available")
            elif(status_code) ==6: ## Data is corrupted or no longer available (cache cleaned, etc)
                print("job running for too long time")
                failed_case.append("job running for too long time")
        except:
            print("Something else went wrong in CDD query")
                        
        
    with open("CDD_domain_additional_bp.txt", "w") as fh:
        for line in CDD_report:
            for i in (range(0,len(line))):
                if(re.search(r'Q#', line[i])): ## chekc if the output is correct, successful query start with Q#
                    fh.write(str(line[i])+"\n")  
#                        print(line[i])
        
    with open("CDD_domain_failed_case_additional_bp.txt", "w") as fh:
        for line in failed_case:
            fh.write(line[0]+"\n")    
    

##
##Search PDB bank with text search using keywor: Histone, cenp-a,H3,H4,H2a,H2b,H1,H5
##retrieve the PDB IDs of the structures related to histones.
first_layer = False
if __name__ == "__main__":
    if (first_layer == False):
        pdb_list = list()
        chains_report =list()
        EntityIDs_list = list()
        chains_report.append("structureId\tchainId\tcompound\tsource\tuniprotAcc\tuniprotRecommendedName\tentityMacromoleculeType")
        search_text("histone",pdb_list)
        search_text("cenp",pdb_list)    
        search_text("H1",pdb_list)    
        search_text("H2a",pdb_list)    
        search_text("H2b",pdb_list)
        search_text("H3",pdb_list) 
        search_text("H4",pdb_list) 
        search_text("H5",pdb_list) 
        search_text("nucleosome",pdb_list) 
    
## remove duplicates and sort the PDB IDs
        pdb_list=list(set(pdb_list)) # remove duplicates 
        pdb_list.sort() 
        print(len(pdb_list))
        with open("PDB_list.txt", "w") as fh:
            for pdb_id in pdb_list:
                try:
                    get_entityIDs(pdb_id,EntityIDs_list)
                    fh.write(pdb_id+"\n") 
                    print(len(EntityIDs_list))
                except:
                    print("Something else went wrong with get_entityIDs")
                    continue
                   
        for EntityIDs in EntityIDs_list:
            try:
                generate_customReport(EntityIDs,chains_report)  
                print(len(chains_report)) 
            except:
                print("Something else went wrong with generate_customReport")
                continue
        
        with open("chains.csv", "w") as fh:
            chains_report = list(filter(None, chains_report))
            for chain in chains_report:
                fh.write(chain+"\n") 
        get_CDD_domain_hit("chains.csv")  
        
        
        
    elif (first_layer == True):
        pdb_list = list()
        chains_report =list()
        EntityIDs_list = list()
        chains_report.append("structureId\tchainId\tcompound\tsource\tuniprotAcc\tuniprotRecommendedName\tentityMacromoleculeType")        
        
        with open("uniprots_bp.txt", "r") as fh:
            for line in fh:
                uniprot_ID = line.strip().upper()
                search_uniprot(uniprot_ID,pdb_list)
             
## remove duplicates and sort the PDB IDs
        pdb_list=list(set(pdb_list)) # remove duplicates 
        pdb_list.sort() 
#        print(len(pdb_list))
        with open("PDB_list.txt", "w") as fh:
            for pdb_id in pdb_list:
                get_entityIDs(pdb_id,EntityIDs_list)
                fh.write(pdb_id+"\n") 
                print(len(EntityIDs_list))
                               
        for EntityIDs in EntityIDs_list:
            generate_customReport(EntityIDs,chains_report)  
            print(len(chains_report)) 
        
        with open("chains.csv", "w") as fh:
            chains_report = list(filter(None, chains_report))
            for chain in chains_report:
                fh.write(chain+"\n") 
        get_CDD_domain_hit("chains.csv")          

    
        










              
                
                
