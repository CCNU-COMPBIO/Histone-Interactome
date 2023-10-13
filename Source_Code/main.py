#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 12:02:07 2019

@author: pengy10
"""

import re
import csv
import Query_PDB_new as Query_PDB
import os

#PATH = "./"
#PATH = "/Users/pengy10/Desktop/Interactome/Interactome_script/new_scripts/all_histone_DNA_interactions_mapping_PDB_numbering/Interfaces/"
PATH = "/media/dell/DiskF/interactomes/pipeline/Interactome/Workflow/Interfaces/"
CHAIN_FILE = "./chains_additional_bp.csv"
#CHAIN_FILE = "tempChains.csv"
PDB_LIST = "./PDB_list_additional_bp.txt"
#PDB_LIST = "tempID.txt"


#PARAMETERS:
#fils is a string with path to the file to be checked

#RESULTS:
#Returns 1 or 0 depending on file existence 


def file_check(file):
    
    try:
        open(file, "r")
        return 1
    
    except IOError:
        print("Error: " + file + " does not appear to exist.")
        return 0

#%%

#PARAMETERS:
#name is a string with the name of chain
#typeCount is a 2-element list with the first one being an empty string and the second - 0 (iteger)

#RESULTS:
#Updates the typeCount array: the first element = general histone type(s); the second element = (should be) number of histones in chain

def is_histone(name, typeCount):  ## for uniprotRecommendedName
    typeCount[1] = 0 ## count number of histone types in a single chain 

    if(re.search(r'histone', name, re.I)):
       
        if(not re.search(r'chaperone|ase|ing|hsh49|rep|alpha|thioredoxin|bombinin|chromosomal|nucleoprotein|envelope|snrnp', name, re.I)):
            
#            typeCount[1] = 1 #adds the number of histones in chain  (should be changed to actual number of histones in chain!!!)

            if(re.search(r'h2a|histone 2a', name, re.I)):
                typeCount[0] += 'H2A||'
                typeCount[1] = typeCount[1] +1
                
            elif(re.search(r'h2b|histone 2b', name, re.I)):
                typeCount[0] += 'H2B||'
                typeCount[1] = typeCount[1] +1
                    
            elif(re.search(r'h3|histone 3', name, re.I)):
                typeCount[0] += 'H3||'
                typeCount[1] = typeCount[1] +1
                
            elif(re.search(r'h4|histone 4', name, re.I)):
                typeCount[0] += 'H4||'
                typeCount[1] = typeCount[1] +1
                    
            elif(re.search(r'h1|h5|histone 1|histone 5', name, re.I)):
                typeCount[0] += 'H1||'  
                typeCount[1] = typeCount[1] +1
                    
#            elif(re.search(r'h5', name, re.I)):
#                typeCount[0] += 'H5|'
#                typeCount[1] = typeCount[1] +1

#            elif(re.search(r'arch', name, re.I)):
#                typeCount[0] += 'archaeal histone|'
#                typeCount[1] = typeCount[1] +1
                
            else:
                typeCount[0] = 'some histone||'
                typeCount[1] = typeCount[1] +1


#%%

#PARAMETERS:
#name is a string with the name of chain
#typeCount is a 2-element list with the first one being an empty string and the second - 0 (iteger)

#RESULTS:
#Updates the typeCount array: the first element = general histone type(s); the second element = (should be) number of histones in chain

def is_histone2(name, typeCount): ## for compound name
    typeCount[1] = 0
    if(re.search(r'histone|h3|h4|h2a|h2b|lys\(ac\)', name, re.I)):       
        if(not re.search(r'chaperone|ase|ing|fasciclin|rna|chain|region|fold|llama|anti|inhibitor|nanobody|nucleoprotein|regulatory|insulin|nr1h|upf0052|hiv|affitin|envelope|ddef|bh3|hsh49|1h4i|h2afy|h2a1', name, re.I)):
            
#            typeCount[1] = 1 #adds the number of histones in chain  (should be changed to actual number of histones in chain!!!)

            if(re.search(r'h2a|histone 2a', name, re.I)):
                typeCount[0] += 'H2A||'
                typeCount[1] = typeCount[1] +1
                
            elif(re.search(r'h2b|histone 2b', name, re.I)):
                typeCount[0] += 'H2B||'
                typeCount[1] = typeCount[1] +1
                    
            elif(re.search(r'h4|histone 4', name, re.I)):
                typeCount[0] += 'H4||'
                typeCount[1] = typeCount[1] +1
                    
            elif(re.search(r'h1|h5|histone 1|histone 5', name, re.I)):
                typeCount[0] += 'H1||'  
                typeCount[1] = typeCount[1] +1
                    
#            elif(re.search(r'h5', name, re.I)):
#                typeCount[0] += 'H5|'
#                typeCount[1] = typeCount[1] +1

#            elif(re.search(r'arch', name, re.I)):
#                typeCount[0] += 'archaeal histone|'
#                typeCount[1] = typeCount[1] +1
                
            elif(re.search(r'h3|histone 3', name, re.I)):
                typeCount[0] += 'H3||'  
                typeCount[1] = typeCount[1] +1
                
            else:
                typeCount[0] = 'some histone||'
                typeCount[1] = typeCount[1] +1


#%%

#PARAMETERS: 
#pdbList is a text file with a header and one column PDB
#files is a list
#parameter is a string, either 'mapping' or 'interface' depending on desired results

#RESULTS:
#A list of absolute paths to either mapping files or interface files as it is stored on local NCBI machines


def get_files(pdbList, files, parameter):
    
    with open(pdbList, 'r') as pfh:
        
        if(parameter == 'mapping'):
            
            for line in pfh:
                line = line.strip().lower()
                folder = line[1] + line[2] 
                files.append(PATH + folder + '/' + line + '_chain_protein_mapping.tab')
                
        elif(parameter == 'interface'): 
            
            for line in pfh:
                line = line.strip().lower()
                folder = line[1] + line[2]
                files.append(PATH + folder + '/' + line + '_atomic_contacts_5.0A.tab')


#%%

#PARAMETERS: 
#pdb is a string with pdb id
#parameter is a string, either 'mapping' or 'interface' depending on desired results

#RESULTS:
#An absolute paths to either mapping file or interface file as it is stored on local NCBI machines


def get_file(pdb, parameter):
    
    if(parameter == 'mapping'):
        folder = pdb[1] + pdb[2] 
        file = (PATH + folder + '/' + pdb + '_chain_protein_mapping.tab')
        return file  
        
    elif(parameter == 'interface'): 
        folder = pdb[1] + pdb[2]
        file = (PATH + folder + '/' + pdb + '_atomic_contacts_5.0A.tab')
        return file


#%%

#PARAMETERS:
#cFile is tab-separated file with a header and 4 columns: pdb, chain, uniprot, name
#dictionary is nested with the innermost dict being dictionary['pdb'] = {}

#RESULTS: 
#Example: {1alq : {'G': 'A|p84233|Histone H3.2|H3|Xenopus laevis|nucleosome:1|bp:1|1|1|135'}}

def get_chain_dictionaries(cFile, dictionary): 
    
    #with open(cFile, 'r') as cfh:
        #cfh.readline()

    histoneCount = {} #is used to count number of histones in a structure!!!!!!!

    tempDict = {}
    tempDict['pdb'] = {}

    tempDict2 = {}
    tempDict2['pdb'] = {}
    
    mappingFiles = [] 
    get_files(PDB_LIST, mappingFiles, 'mapping')


    for file in mappingFiles:

        try: #adds a pdb entry to the dict only if mapping file exists

            with open(file, 'r') as mfh:
                mfh.readline() #skips header
                pdb = file.split('/')[-1].split('_', 1)[0]

                for mLine in mfh:
                    chainPair = mLine.split('\t')
                    alexChain = chainPair[0]
                    myChain = chainPair[1]
                    alignment1 = int(chainPair[7]) - int(chainPair[3]) # align to uniprot sequence number
                    alignment2 = int(chainPair[9]) - int(chainPair[3]) # align to PDB sequence number
                    mmCIFstart = int(chainPair[3])
                    mmCIFend = int(chainPair[5])
                    
                    if(pdb in tempDict):
                        if(alexChain in tempDict[pdb]):
                            tempDict2[pdb][alexChain].extend([alignment1,alignment2, mmCIFstart, mmCIFend])
                        else:
                            tempDict[pdb][alexChain] = myChain
                            tempDict2[pdb][alexChain] = [alignment1,alignment2, mmCIFstart, mmCIFend]
                        
                    else:
                        tempDict[pdb] = {alexChain : myChain}
                        tempDict2[pdb] = {alexChain : [alignment1,alignment2, mmCIFstart, mmCIFend]}

        except IOError:
            pass
            #print("Error: " + mappingFile + " does not appear to exist.")

    chains_file = csv.DictReader(open(cFile), delimiter='\t') 

    for cLine in chains_file:

        pdb = cLine["structureId"].lower()
        if(pdb in tempDict): #continues only if a mapping file exists
            chain = cLine["chainId"]
            organism = cLine["source"]
            uniprot = cLine["uniprotAcc"].lower()
            name = cLine["uniprotRecommendedName"]
            chaintype = cLine["entityMacromoleculeType"]

            histoneTypeAndCount = ['', 0]

            if(name!="NA"):
                is_histone(name, histoneTypeAndCount) #checks whether the name looks like a histone!!!!!!!!!

            else:
                name = cLine["compound"]
                is_histone2(name, histoneTypeAndCount)

            tempType = histoneTypeAndCount[0]
            tempCount = histoneTypeAndCount[1]

            #######################
            if(tempCount): #if the chain is a [part of a] histone

                if(pdb in histoneCount):

                    if(tempType not in histoneCount[pdb]):
                        histoneCount[pdb].append(tempType) #!!!!!!

                else:
                    histoneCount[pdb] = [tempType]       

            try: #adds a chain entry to the dict only if there exists a corresponding chain in the mapping file
                alexChain = list(tempDict[pdb].keys())[list(tempDict[pdb].values()).index(chain)]

                if(pdb in dictionary):

                    if((tempCount==1) and (re.search(r'Protein', chaintype, re.I))): #checks if chain is a histone!!!!c and exclude the non-protein chain
                        dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '||' + uniprot + '||' + name + '||' + tempType + organism + '||'  #!!!!
                    elif((tempCount>1) and (re.search(r'Protein', chaintype, re.I))): ## single chain classified as mutiple histone type (manully check required)
                        dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '||' + uniprot + '||' + name + '||' + 'mutiple_histone_type||' + organism + '||' #!!!!

                    elif(re.search(r'Protein', chaintype, re.I)): #!!!! exclude the non-protein chain
                        dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '||' + uniprot + '||' + name + '||' + 'other||' + organism  + '||' #!!!!
                    
                    elif(re.search(r'DNA', chaintype, re.I)): #!!!! exclude the non-protein chain
                        dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '||' + uniprot + '||' + name + '||' + 'DNA||' + organism  + '||' #!!!!

                else:

                    if((tempCount==1) and (re.search(r'Protein', chaintype, re.I))): #checks if chain is a histone!!!!
                        dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '||' + uniprot + '||' + name + '||' + tempType + organism + '||'} #!!!!
                    elif((tempCount>1) and (re.search(r'Protein', chaintype, re.I))):
                        dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '||' + uniprot + '||' + name + '||' + 'mutiple_histone_type||' + organism + '||'} #!!!!
                    elif(re.search(r'Protein', chaintype, re.I)): #!!!!
                        dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '||' + uniprot + '||' + name + '||' + 'other||' + organism + '||'} #!!!!
                    elif(re.search(r'DNA', chaintype, re.I)): #!!!!
                        dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '||' + uniprot + '||' + name + '||' + 'DNA||' + organism + '||'} #!!!!


            except ValueError:
                #print("Error: " + str(ValueError) + ", in " + pdb)               
                pass


    with open('temp.tsv', 'w') as rfh:
        for structure in list(dictionary): 

            if(structure in histoneCount):
                uniqueHistoneNum = len(histoneCount[structure])
                partnerFlag = 0

                if(uniqueHistoneNum > 2): #checks if pdb has at least 3 different histones ~ is a nucleosome!!!!!!!

                    for chain in dictionary[structure]:
                        chainFields = dictionary[structure][chain].split('||')

                        chainType = chainFields[3]  ##

                        dictionary[structure][chain] += 'nucleosome:1||' #!!!!!!

                        if(partnerFlag == 0 and chainType == 'other'):
                            partnerFlag = 1
                            
                        if(partnerFlag == 0 and chainType == 'DNA'):
                            partnerFlag = 1                           
                        

                    if(partnerFlag == 0):
                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:0||'
                            for element in tempDict2[structure][chain]:
                                dictionary[structure][chain] += str(element)
                                dictionary[structure][chain] += '||'
                            dictionary[structure][chain] += str(len(tempDict2[structure][chain]) / 4) ## single chain could have mutiply mapping, we need to check this.
                            rfh.write(structure + "||" + str(dictionary[structure][chain]) + '\n')
                            
                            #rfh.write(structure + '\t' + 'nucleosome' + '\t' + 'no' + '\n')
                    else:
                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:1||'
                            for element in tempDict2[structure][chain]:
                                dictionary[structure][chain] += str(element)
                                dictionary[structure][chain] += '||'
                            dictionary[structure][chain] += str(len(tempDict2[structure][chain]) / 4)
                            rfh.write(structure + "||" +  str(dictionary[structure][chain]) + '\n')
                            #rfh.write(structure + '\t' + 'nucleosome' + '\t' + 'yes' + '\n')
                else: #!!!!!!

                    for chain in dictionary[structure]: #!!!!!
                        chainFields = dictionary[structure][chain].split('||')

                        chainType = chainFields[3]

                        dictionary[structure][chain] += 'nucleosome:0||' #!!!!!! 

                        if(partnerFlag == 0 and chainType == 'other'):
                            partnerFlag = 1
                            
                        if(partnerFlag == 0 and chainType == 'DNA'):
                            partnerFlag = 1

                    if(partnerFlag == 0):
                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:0||'
                            for element in tempDict2[structure][chain]:
                                dictionary[structure][chain] += str(element)
                                dictionary[structure][chain] += '||'
                            dictionary[structure][chain] += str(len(tempDict2[structure][chain]) / 4)
                            rfh.write(structure + "||" +  str(dictionary[structure][chain]) + '\n')
                            #rfh.write(structure + '\t' + 'histone' + '\t' + 'no' + '\n')
                    else:
                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:1||'
                            for element in tempDict2[structure][chain]:
                                dictionary[structure][chain] += str(element)
                                dictionary[structure][chain] += '||'
                            dictionary[structure][chain] += str(len(tempDict2[structure][chain]) / 4)
                            rfh.write(structure + "||" +  str(dictionary[structure][chain]) + '\n')
                            #rfh.write(structure + '\t' + 'histone' + '\t' + 'yes' + '\n')
            else:
                for chain in dictionary[structure]:
                    chainFields = dictionary[structure][chain].split('||')

                    chainType = chainFields[3]  ##

                    dictionary[structure][chain] += 'other:0||' #!!!!!!
                            
                    #rfh.write(structure + '\t' + 'nucleosome' + '\t' + 'no' + '\n')

                for chain in dictionary[structure]:
                    dictionary[structure][chain] += 'bp:1||'
                    for element in tempDict2[structure][chain]:
                        dictionary[structure][chain] += str(element)
                        dictionary[structure][chain] += '||'
                    dictionary[structure][chain] += str(len(tempDict2[structure][chain]) / 4)
                    rfh.write(structure + "||" +  str(dictionary[structure][chain]) + '\n')
                #del dictionary[structure]


#%%

#PARAMETERS:
#interfaceFiles is a list of strings containing names of interface files
#chainDictionary is a dictionary produced by 'get_chain_dictionaries'
#interfaceDictonary is a nested dictionary with the innermost dict being interfaceDictionary['uniprotPair'] = {}, where the values are lists of the next form [chain pair data, residue number, pdb IDs] 
#Example: 'P84233@P62799': {'44': ['A|P84233|Histone H3.2|h3|1@B|P62799|Histone H4|h4|1$E|P84233|Histone H3.2|h3|1@F|P62799|Histone H4|h4|1', 17, '1zla$q5cl']
#note that data fields of a chain is separated by |, @ separates chain from binding partner chain, $ separates chain pairs and pdb structures
#A|P84233|Histone H3.2|H3|Xenopus laevis|nucleosome:1|bp:1

## Falg =0: find the structures for histones with binding partners
## Flag =1: find the fistr layer of binding parterns 
def residue_count(interfaceFiles, chainDictionary, interfaceDictionary,flag_first_layer = False ): 
    with open('results_all_histone_PPIs_additional_bp.tsv', 'w') as rfh:    
        for file in interfaceFiles:
            pdb = file.split('/')[-1].split('_', 1)[0] 

            try:
            
                with open (file, 'r') as ifh:
                    ifh.readline()
                    interface_residue_pairs = list()
                    uniq_list = list()
                    interface_residue_pairs = [[line.split('\t')[0],line.split('\t')[1],line.split('\t')[2],line.split('\t')[4],line.split('\t')[5],line.split('\t')[6]]for line in ifh]
                    if(len(interface_residue_pairs) <= 500000): ##to exclude large structures like HIV virus and save time
                        for sublist in interface_residue_pairs:
                            if sublist not in uniq_list:
                               uniq_list.append(sublist)

#                with open('results_all.tsv', 'a') as rfh:
                    
                    for lineFields in uniq_list:

                        chain1 = lineFields[0].split('_', 1)[0] # the split part treats biological assembly chains as separate chains ???
                        chain2 = lineFields[3].split('_', 1)[0]
                        
                        resnum1 = lineFields[2]    ##interfacial residue position and residue name in Alex's atomic contact file 
                        resnum2 = lineFields[5]
                        resname1 = lineFields[1]
                        resname2 = lineFields[4]
                        
                        if (chain1 in chainDictionary[pdb] and chain2 in chainDictionary[pdb]):
                            fields1 = chainDictionary[pdb][chain1].split('||')
                            fields2 = chainDictionary[pdb][chain2].split('||')
                        else:
                            continue

                        myChain1 = fields1[0]
                        myChain2 = fields2[0]
                        
                        type1 = fields1[3]
                        type2 = fields2[3]
                        
                        flag_mapping1 = int(float(fields1[-1]))  ## Flag for mapping file
                        flag_mapping2 = int(float(fields2[-1]))  ## Flag for mapping file
                        

                        
                        if (flag_mapping1 == 1 and flag_mapping2 ==1): ## Single chain have single mapping in the Alex's mapping file 
                        #nucleosome = int(fields1[5].split(':')[1])
#                        bp = int(fields1[6].split(':')[1])
                            align1_uniprot = int(fields1[7])
                            align2_uniprot = int(fields2[7])
                            align1_pdb = int(fields1[8])
                            align2_pdb = int(fields2[8])
                        
                            mmCIFstart1 = int(fields1[9])
                            mmCIFstart2 = int(fields2[9])
                        
                            mmCIFend1 = int(fields1[10])
                            mmCIFend2 = int(fields2[10])
                        
                            uniprotStart1 = mmCIFstart1 + align1_uniprot
                            uniprotStart2 = mmCIFstart2 + align2_uniprot
                        
                            uniprotEnd1 = mmCIFend1 + align1_uniprot
                            uniprotEnd2 = mmCIFend2 + align2_uniprot
                        elif (flag_mapping1 > 1 and flag_mapping2 ==1): ## Single chain have mutiple mapping in the Alex's mapping file 
                            for i in range(1,flag_mapping1 + 1):
                                if (int(resnum1) >= int(fields1[9+(i-1)*4])) and int(resnum1) <= int(fields1[10+(i-1)*4]): ## select the mapping chains where the residue corresponds to 
                                    
                                    align1_uniprot = int(fields1[7+(i-1)*4])
                                    align2_uniprot = int(fields2[7])
                                    align1_pdb = int(fields1[8+(i-1)*4])
                                    align2_pdb = int(fields2[8])
                                    
                                    mmCIFstart1 = int(fields1[9+(i-1)*4])
                                    mmCIFstart2 = int(fields2[9])
                        
                                    mmCIFend1 = int(fields1[10+(i-1)*4])
                                    mmCIFend2 = int(fields2[10])
                        
                                    uniprotStart1 = mmCIFstart1 + align1_uniprot
                                    uniprotStart2 = mmCIFstart2 + align2_uniprot
                        
                                    uniprotEnd1 = mmCIFend1 + align1_uniprot
                                    uniprotEnd2 = mmCIFend2 + align2_uniprot 
                                    
                        elif (flag_mapping2 > 1 and flag_mapping1 ==1): ## Single chain have mutiple mapping in the Alex's mapping file 
                            for i in range(1,flag_mapping2 + 1):
                                if (int(resnum2) >= int(fields2[9+(i-1)*4])) and int(resnum2) <= int(fields2[10+(i-1)*4]):
                                    
                                    align1_uniprot = int(fields1[7])
                                    align2_uniprot = int(fields2[7+(i-1)*4])
                                    align1_pdb = int(fields1[8])
                                    align2_pdb = int(fields2[8+(i-1)*4])
                                    
                                    mmCIFstart1 = int(fields1[9])
                                    mmCIFstart2 = int(fields2[9+(i-1)*4])
                        
                                    mmCIFend1 = int(fields1[10])
                                    mmCIFend2 = int(fields2[10+(i-1)*4])
                        
                                    uniprotStart1 = mmCIFstart1 + align1_uniprot
                                    uniprotStart2 = mmCIFstart2 + align2_uniprot
                        
                                    uniprotEnd1 = mmCIFend1 + align1_uniprot
                                    uniprotEnd2 = mmCIFend2 + align2_uniprot 
                                
                        elif (flag_mapping1 > 1 and flag_mapping2 > 1): ## Single chain have mutiple mapping in the Alex's mapping file 
                            for i in range(1,flag_mapping1 + 1):
                                if (int(resnum1) >= int(fields1[9+(i-1)*4])) and (int(resnum1) <= int(fields1[10+(i-1)*4])):
                                    align1_uniprot = int(fields1[7+(i-1)*4])
                                    align1_pdb = int(fields1[8+(i-1)*4])
                                    
                                    mmCIFstart1 = int(fields1[9+(i-1)*4])
                                    mmCIFend1 = int(fields1[10+(i-1)*4])
                                    uniprotStart1 = mmCIFstart1 + align1_uniprot
                                    uniprotEnd1 = mmCIFend1 + align1_uniprot
                                    
                            for i in range(1,flag_mapping2 + 1):
                                if (int(resnum2) >= int(fields2[9+(i-1)*4])) and (int(resnum2) <= int(fields2[10+(i-1)*4])):
                                    align2_uniprot = int(fields2[7+(i-1)*4])
                                    align2_pdb = int(fields2[8+(i-1)*4])
                                    
                                    mmCIFstart2 = int(fields2[9+(i-1)*4])
                                    mmCIFend2 = int(fields2[10+(i-1)*4])
                                    uniprotStart2 = mmCIFstart2 + align2_uniprot
                                    uniprotEnd2 = mmCIFend2 + align2_uniprot 
                                    
                        ###### output the information of differnt type of histone interactions
                        if(((type1 != 'other' and type1 != 'DNA') or (type1=='DNA' and type2=='other')) and flag_first_layer == False):
                        #if(type1 != 'other'):
                            #if((type1 == 'H1' and 1 <= uniprotStart1 <= 229 and 1 <= uniprotEnd1 <= 229) or (type1 == 'H2A' and 1 <= uniprotStart1 <= 385 and 1 <= uniprotEnd1 <= 385) or (type1 == 'H2B' and 1 <= uniprotStart1 <= 133 and 1 <= uniprotEnd1 <= 133) or (type1 == 'H3' and 1 <= uniprotStart1 <= 184 and 1 <= uniprotEnd1 <= 184) or (type1 == 'H4' and 1 <= uniprotStart1 <= 103 and 1 <= uniprotEnd1 <= 103)):
                                #if(int(resnum1) >= mmCIFstart1 and int(resnum1) <= mmCIFend1):
                            with open('./CDD_domain.txt', 'r') as hfh:
                                hfh.readline()
                                domainFlag = 0
                                
                                for line in hfh:
                                    fields = line.split('\t')


                                    pdb2 = fields[0].split(' - ')[1][0:4].lower() ###
                                    chain = fields[0].split(' - ')[1][5]  ###


                                    start = int(fields[3])
                                    end = int(fields[4])
                                    domtype = fields[1]
                                    name = fields[8]

                                    if(pdb2 == pdb and chain == myChain2):

                                        if(domtype == 'specific'):
                                            if(int(resnum2) >= int(start) and int(resnum2) <= int(end)):
                                                hfields = chainDictionary[pdb][chain1].split('||')
                                                pfields = chainDictionary[pdb][chain2].split('||')
                                                
                                                alignedRes_uniprot1 = int(resnum1) + align1_uniprot ## residue number in unprot sequence
                                                alignedRes_pdb1 = int(resnum1) + align1_pdb ## residue number in PDB file
                                                alignedRes_uniprot2 = int(resnum2) + align2_uniprot
                                                alignedRes_pdb2 = int(resnum2) + align2_pdb
                                                
                                                domainFlag = 1
                                                rfh.write(hfields[1] +'\t' + hfields[0] + "\t" + hfields[2] + '\t' + hfields[4] + '\t' + hfields[3] + '\t' + str(alignedRes_uniprot1) + '\t' + str(alignedRes_pdb1) + '\t' + resname1 + '\t' + pfields[1] + '\t' + pfields[0] + '\t' + pfields[2] + '\t' + name + '\t' + pfields[4] + '\t' + str(alignedRes_uniprot2) + '\t' + str(alignedRes_pdb2) + '\t' + resname2 + '\t' + pdb + '\t' + hfields[5] + '\t' + pfields[3] + '\n')
                                                #print(hfields[2] + '\t' + hfields[3] + '\t' + hfields[4] + '\t' + resnum1 + '\t' + pfields[1] + '\t' + pfields[2] + '\t' + pfields[4] + '\t' + name)
                                hfh.seek(0)
                                if(not domainFlag):  ### domainFlag == 0
                                    hfields = chainDictionary[pdb][chain1].split('||')
                                    pfields = chainDictionary[pdb][chain2].split('||')
                                    alignedRes_uniprot1 = int(resnum1) + align1_uniprot ## residue number in unprot sequence
                                    alignedRes_pdb1 = int(resnum1) + align1_pdb ## residue number in PDB file
                                    alignedRes_uniprot2 = int(resnum2) + align2_uniprot
                                    alignedRes_pdb2 = int(resnum2) + align2_pdb
                                    rfh.write(hfields[1] +'\t' + hfields[0] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[3] + '\t' + str(alignedRes_uniprot1) + '\t' + str(alignedRes_pdb1) + '\t' + resname1 + '\t' + pfields[1] + '\t' + pfields[0] + '\t' + pfields[2] + '\t' + 'NA' + '\t' + pfields[4] + '\t' + str(alignedRes_uniprot2) + '\t' + str(alignedRes_pdb2)  +'\t' + resname2 + '\t' + pdb + '\t' + hfields[5] + '\t' + pfields[3] + '\n')
             

                        if(((type2 != 'other' and type2 != 'DNA') or (type2=='DNA' and type1=='other')) and  flag_first_layer == False):
                        #if(type2 != 'other' and type2 != 'DNA'):
                        #if(type2 != 'other'):                            
                            #if((type2 == 'H1' and 1 <= uniprotStart2 <= 229 and 1 <= uniprotEnd2 <= 229) or (type2 == 'H2A' and 1 <= uniprotStart2 <= 385 and 1 <= uniprotEnd2 <= 385) or (type2 == 'H2B' and 1 <= uniprotStart2 <= 133 and 1 <= uniprotEnd2 <= 133) or (type2 == 'H3' and 1 <= uniprotStart2 <= 184 and 1 <= uniprotEnd2 <= 184) or (type2 == 'H4' and 1 <= uniprotStart2 <= 103 and 1 <= uniprotEnd2 <= 103)):
                                #if(int(resnum2) >= mmCIFstart2 and int(resnum2) <= mmCIFend2):
                            with open('./CDD_domain.txt', 'r') as hfh:
                                hfh.readline()
                                domainFlag = 0
                                for line in hfh:
                                    fields = line.split('\t')


                                    pdb2 = fields[0].split(' - ')[1][0:4].lower() ###
                                    chain = fields[0].split(' - ')[1][5]  ###

                                    start = int(fields[3])
                                    end = int(fields[4])
                                    domtype = fields[1]
                                    name = fields[8]
                                    if(pdb2 == pdb and chain == myChain1):

                                        if(domtype == 'specific') :
                                            if(int(resnum1) >= int(start) and int(resnum1) <= int(end)):
                                                hfields = chainDictionary[pdb][chain2].split('||')
                                                pfields = chainDictionary[pdb][chain1].split('||')
                                                alignedRes_uniprot2 = int(resnum2) + align2_uniprot
                                                alignedRes_pdb2 = int(resnum2) + align2_pdb
                                                alignedRes_uniprot1 = int(resnum1) + align1_uniprot ## residue number in unprot sequence
                                                alignedRes_pdb1 = int(resnum1) + align1_pdb ## residue number in PDB file
                                                
                                                domainFlag = 1
                                                rfh.write(hfields[1] +'\t' + hfields[0] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[3] + '\t' + str(alignedRes_uniprot2) + '\t' + str(alignedRes_pdb2) + '\t' + resname2 + '\t' + pfields[1] + '\t' + pfields[0] + '\t' + pfields[2] + '\t' + name + '\t' + pfields[4] + '\t' + str(alignedRes_uniprot1) + '\t' + str(alignedRes_pdb1)  + '\t' + resname1 + '\t' + pdb + '\t' + hfields[5] + '\t' + pfields[3] + '\n')

                                                #print(hfields[2] + '\t' + hfields[3] + '\t' + hfields[4] + '\t' + resnum2 + '\t' + pfields[1] + '\t' + pfields[2] + '\t' + pfields[4] + '\t' + name)
                                hfh.seek(0)
                                if(not domainFlag):
                                    hfields = chainDictionary[pdb][chain2].split('||')
                                    alignedRes_uniprot2 = int(resnum2) + align2_uniprot
                                    alignedRes_pdb2 = int(resnum2) + align2_pdb
                                    alignedRes_uniprot1 = int(resnum1) + align1_uniprot ## residue number in unprot sequence
                                    alignedRes_pdb1 = int(resnum1) + align1_pdb ## residue number in PDB file
                                    
                                    pfields = chainDictionary[pdb][chain1].split('||')
                                    rfh.write(hfields[1] +'\t' + hfields[0] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[3] + '\t' + str(alignedRes_uniprot2) + '\t' + str(alignedRes_pdb2) + '\t' + resname2 + '\t' + pfields[1] + '\t' + pfields[0] + '\t' + pfields[2] + '\t' + 'NA' + '\t' + pfields[4] + '\t' + str(alignedRes_uniprot1) + '\t' + str(alignedRes_pdb1) + '\t' + resname1 + '\t' + pdb + '\t' + hfields[5] + '\t' + pfields[3] + '\n')
                        
                        if(type1 == 'other' and type2 == 'other' and flag_first_layer == True):
                            with open('./CDD_domain.txt', 'r') as hfh:
                                hfh.readline()
                                domainFlag1 = 0
                                domainFlag2 = 0
                                for line in hfh:
                                    fields = line.split('\t')
                                    pdb2 = fields[0].split(' - ')[1][0:4].lower() ###
                                    chain = fields[0].split(' - ')[1][5]  ###
                                    start = int(fields[3])
                                    end = int(fields[4])
                                    domtype = fields[1]
                                    name = fields[8]
                                    if(pdb2 == pdb and chain == myChain1):
                                        if(domtype == 'specific') :
                                            if(int(resnum1) >= int(start) and int(resnum1) <= int(end)):
                                                hfields = chainDictionary[pdb][chain2].split('||')
                                                pfields = chainDictionary[pdb][chain1].split('||')
                                                alignedRes_uniprot2 = int(resnum2) + align2_uniprot
                                                alignedRes_pdb2 = int(resnum2) + align2_pdb
                                                alignedRes_uniprot1 = int(resnum1) + align1_uniprot ## residue number in unprot sequence
                                                alignedRes_pdb1 = int(resnum1) + align1_pdb ## residue number in PDB file
                                                domainFlag1 = 1
                                                rfh.write(hfields[1] +'\t' + hfields[0] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[1] + '\t' + str(alignedRes_uniprot2) + '\t' + str(alignedRes_pdb2) + '\t' + resname2 + '\t' + pfields[1]  + '\t' + pfields[0] + '\t' + pfields[2] + '\t' + name + '\t' + pfields[4] + '\t' + str(alignedRes_uniprot1) + '\t' + str(alignedRes_pdb1) + '\t' + resname1 + '\t' + pdb + '\t' + hfields[5] + '\t' + pfields[3] +"_his_bp" + '\n')
                                    elif(pdb2 == pdb and chain == myChain2):
                                        if(domtype == 'specific') :
                                            if(int(resnum2) >= int(start) and int(resnum2) <= int(end)):
                                                hfields = chainDictionary[pdb][chain1].split('||')
                                                pfields = chainDictionary[pdb][chain2].split('||')
                                                alignedRes_uniprot1 = int(resnum1) + align1_uniprot
                                                alignedRes_pdb1 = int(resnum1) + align1_pdb
                                                alignedRes_uniprot2 = int(resnum2) + align2_uniprot
                                                alignedRes_pdb2 = int(resnum2) + align2_pdb
                                                domainFlag2 = 1
                                                rfh.write(hfields[1] +'\t' + hfields[0] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[1] + '\t' + str(alignedRes_uniprot1) + '\t' + str(alignedRes_pdb1) + '\t' + resname1 + '\t' + pfields[1]  + '\t' + pfields[0] + '\t' + pfields[2] + '\t' + name + '\t' + pfields[4] + '\t' + str(alignedRes_uniprot2) + '\t' + str(alignedRes_pdb2) + '\t' + resname2 + '\t' + pdb + '\t' + hfields[5] + '\t' + pfields[3] +"_his_bp" + '\n')

                                if(domainFlag1 ==0): 
                                    hfields = chainDictionary[pdb][chain2].split('||')
                                    pfields = chainDictionary[pdb][chain1].split('||')
                                    alignedRes_uniprot2 = int(resnum2) + align2_uniprot
                                    alignedRes_pdb2 = int(resnum2) + align2_pdb
                                    alignedRes_uniprot1 = int(resnum1) + align1_uniprot
                                    alignedRes_pdb1 = int(resnum1) + align1_pdb
                                    rfh.write(hfields[1] +'\t' + hfields[0] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[1] + '\t' + str(alignedRes_uniprot2) + '\t' + str(alignedRes_pdb2) + '\t' + resname2 + '\t' + pfields[1]  + '\t' + pfields[0] + '\t' + pfields[2] + '\t' + 'NA' + '\t' + pfields[4] + '\t' + str(alignedRes_uniprot1) + '\t' + str(alignedRes_pdb1) + '\t' + resname1 + '\t' + pdb + '\t' + hfields[5] + '\t' + pfields[3] +"_his_bp" + '\n')
                            
                                if(domainFlag2 ==0):
                                    hfields = chainDictionary[pdb][chain1].split('||')
                                    pfields = chainDictionary[pdb][chain2].split('||')
                                    alignedRes_uniprot1 = int(resnum1) + align1_uniprot
                                    alignedRes_pdb1 = int(resnum1) + align1_pdb
                                    alignedRes_uniprot2 = int(resnum2) + align2_uniprot
                                    alignedRes_pdb2 = int(resnum2) + align2_pdb
                                    rfh.write(hfields[1] +'\t' + hfields[0] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[1] + '\t' + str(alignedRes_uniprot1) + '\t' + str(alignedRes_pdb1) + '\t' + resname1 + '\t' + pfields[1]  + '\t' + pfields[0] + '\t' + pfields[2] + '\t' + 'NA' + '\t' + pfields[4] + '\t' + str(alignedRes_uniprot2) + '\t' + str(alignedRes_pdb2) + '\t' + resname2 + '\t' + pdb + '\t' + hfields[5] + '\t' + pfields[3] +"_his_bp" + '\n')
                                
                        
            except (IOError, KeyError) as e:
                #print("Error: " + file + " does not appear to exist.")
                print(e)
                pass


        
#%%
## Output the structures of nucleosome with binding partners; Sine binding parrners can only bind with nucleosomal DNA and these 
## cases are not included in the Alex's protein-protein interaction files. Thus, this function will output the nucleosome structure with binding partner 
## to histone or DNA for generation of interactome of nuclesome
def nucleosome_bp():
    with open('nucleosome_bp_all_PPIs.tsv', 'w') as rfh:
        with open('./temp.tsv', 'r') as hfh:
            hfh.readline()
            for line in hfh:
                fields = line.split('||')
                ptype = fields[4]
                nucleosome_flag = fields[6]
                bp_flag = fields[7]
                if(nucleosome_flag == "nucleosome:1" and bp_flag == "bp:1" and ptype =="other"):
                    rfh.write("nucleosome" + '\t' + fields[2] + '\t' + fields[5] + '\t' + fields[3] + '\t' + fields[0] + '\t' + fields[1] +'\n') 
                
        
    
    
#%%

def main(flag = True):
### search PDB related with histones and stored in pdb_list 

    if (flag == False):
        pdb_list = list()
        chains_report =list()
        EntityIDs_list = list()
        chains_report.append("structureId\tchainId\tcompound\tsource\tuniprotAcc\tuniprotRecommendedName\tentityMacromoleculeType")
        Query_PDB.search_text("histone",pdb_list)
        Query_PDB.search_text("cenp",pdb_list)    
        Query_PDB.earch_text("H1",pdb_list)    
        Query_PDB.search_text("H2a",pdb_list)    
        Query_PDB.search_text("H2b",pdb_list)
        Query_PDB.search_text("H3",pdb_list) 
        Query_PDB.search_text("H4",pdb_list) 
        Query_PDB.search_text("H5",pdb_list) 
        Query_PDB.search_text("nucleosome",pdb_list) 
    
## remove duplicates and sort the PDB IDs
        pdb_list=list(set(pdb_list)) # remove duplicates 
        pdb_list.sort() 
        print(len(pdb_list))
        with open("PDB_list.txt", "w") as fh:
            for pdb_id in pdb_list:
                try:
                    Query_PDB.get_entityIDs(pdb_id,EntityIDs_list)
                    fh.write(pdb_id+"\n") 
                    print(len(EntityIDs_list))
                except:
                    print("Something else went wrong with get_entityIDs")
                    continue
                   
        for EntityIDs in EntityIDs_list:
            try:
                Query_PDB.generate_customReport(EntityIDs,chains_report)  
                print(len(chains_report)) 
            except:
                print("Something else went wrong with generate_customReport")
                continue
        
        with open("chains.csv", "w") as fh:
            chains_report = list(filter(None, chains_report))
            for chain in chains_report:
                fh.write(chain+"\n") 
        Query_PDB.get_CDD_domain_hit("chains.csv") 
        
        
        
    elif (flag == True): ##find the first layer of binding partners.
        """
        pdb_list = list()
        chains_report =list()
        EntityIDs_list = list()
        chains_report.append("structureId\tchainId\tcompound\tsource\tuniprotAcc\tuniprotRecommendedName\tentityMacromoleculeType")        
        
        with open("histone_bp_uniprot_IDs.txt", "r") as fh:
            for line in fh:
                try:                       
                    uniprot_ID = line.strip().upper()
                    Query_PDB.search_uniprot(uniprot_ID,pdb_list)
                except:
                    print("Something else went wrong with search_uniprot")
                    continue
                print(len(pdb_list))
             
## remove duplicates and sort the PDB IDs
        pdb_list=list(set(pdb_list)) # remove duplicates 
        pdb_list.sort() 
        print(pdb_list)
        with open("PDB_list_additional_bp.txt", "w") as fh:
            for pdb_id in pdb_list:
                try:
                    Query_PDB.get_entityIDs(pdb_id,EntityIDs_list)
                    fh.write(pdb_id+"\n") 
                    print(len(EntityIDs_list))
                except:
                    print("Something else went wrong with get_entityIDs")
                    continue
        
        for EntityIDs in EntityIDs_list:
            try:                       
                Query_PDB.generate_customReport(EntityIDs,chains_report)  
                print(len(chains_report)) 
            except:
                print("Something else went wrong with generate_customReport")
                continue
        
        with open("chains_additional_bp.csv", "w") as fh:
            chains_report = list(filter(None, chains_report))
            for chain in chains_report:
                fh.write(chain+"\n") 
        """
        Query_PDB.get_CDD_domain_hit("chains_additional_bp.csv")   
        
    elif (flag == None):
        pass
#    os.system("grep \"Q\#1\" CDD_domain.txt > CDD_domain.tmp")
#    os.system("mv CDD_domain.tmp CDD_domain.txt")
    
            
### generate the histone interactome            
    chainDictionary = {}
    chainDictionary['chain'] = {}
    get_chain_dictionaries(CHAIN_FILE, chainDictionary) ## perform chain mapping with Alex's file 
    nucleosome_bp()  ## get binding paterners interact with nucleosome
    #for pdb in chainDictionary:
     #   for chain in chainDictionary[pdb]:
      #      print(pdb + '\t' + chain + '\t' + str(chainDictionary[pdb][chain]))

    interfaceFiles = []
    get_files(PDB_LIST, interfaceFiles, 'interface')

    interfaceDictionary = {}
    interfaceDictionary['uniprotPair'] = {}
    residue_count(interfaceFiles, chainDictionary, interfaceDictionary,True)
    #for pair in interfaceDictionary:
    #    for a in interfaceDictionary[pair]:
    #        print(pair + '\t' + a + '\t' + str(interfaceDictionary[pair][a]))


#%%

if __name__ == "__main__":
    main(flag = True)

##%
################################################Chain Annotation########################################################
# In the Alex mapping file:  like, 6b3w_chain_protein_mapping.tab
#The chain ID in alex file and PDB is not the same 
# mmCIFstart: Column 4: seq_aln_begin, mmCIFend: Column 6: seq_aln_end
# This is the beginning and ending residue number aligned to the corresponding uniprot in the PDB sequence;
# For example, the pdb sequecne from 10 to 100 can be aligned to unprot, mmCIFstart =10, mmCIFend=100 ( ths is different from the residue number in PDB)

###alignment = db_aln_begin - seq_aln_begin: give the differences between residue number in PDB sequence and number aligned in uniprot sequence
##  db_aln_begin is the residue number aligned in the uniprot sequence, for example, the first residue in PDB sequence is alignment to 10th residue in uniprot sequecen,
## db_aln_begin will be 10

## In the CDD domain files, the domain region number are using the residue number in PDB sequence ( not the unproit residue number or residue number in PDB)
### in the atomic contacts files: like 6b3w_atomic_contacts_5.0A.tab
## the residue number in column 3 and 7 are the number in the PDB sequence, this number is same as the seq_aln_begin in the mapping file and domain region numberiing in the CDD domain files.
########################################################################################################################
