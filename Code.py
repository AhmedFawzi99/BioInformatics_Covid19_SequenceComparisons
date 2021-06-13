#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import csv
import numpy as np


# In[ ]:


#import the data using Bio.Seq Library
seq_file = SeqIO.parse("All Seqs.fasta", "fasta")
seq_file2 = SeqIO.parse("Alignment.clustal_num", "clustal")

seq_count = 0    #Count for the number of Sequences
#Arrays for the Sequence Data
seq_list = []
SequenceName =[]
seq_lengths_list = []

# Arrays for the Neuclotide Base content in each sequence    
AContent=[]
TContent=[]
GContent=[]
CContent=[]
CGContent=[]

# For loop to append the sequence data from the fasta file loaded 
for record in seq_file:
    seq_count+=1
    seq_list.append(record)
    seq_lengths_list.append(len(record.seq))
    


# In[ ]:


# Creating a Dictionary to sort and store the data
Dictionary = { "Sequence Name":[] , "Sequences": [], "A Content" :[], "T Content" :[], "C Content" :[],  
              "G Content":[], 
              "CG Content":[]}


# In[ ]:


# Filling the data in the Dictionary and calculating the percentages required
for i in range(seq_count):
    
    Dictionary["Sequence Name"].append(seq_list[i].name)
    Dictionary["Sequences"].append(seq_list[i].seq)
    Dictionary["A Content"].append(((seq_list[i].seq.count("A"))/len(seq_list[i].seq))*100)
    Dictionary["T Content"].append(((seq_list[i].seq.count("T"))/len(seq_list[i].seq))*100)
    Dictionary["C Content"].append(((seq_list[i].seq.count("C"))/len(seq_list[i].seq))*100)
    Dictionary["G Content"].append(((seq_list[i].seq.count("G"))/len(seq_list[i].seq))*100)
    Dictionary["CG Content"].append(((seq_list[i].seq.count("C")+seq_list[i].seq.count("G"))/len(seq_list[i].seq))*100)



# In[ ]:

#converting the percentages dictionary into a csv file
Data=pd.DataFrame.from_dict(Dictionary)
# Saving the Data as csv 
Data.to_csv('Chemical Constituents.csv')









# In[ ]:


# opening the Allignment file downloaded from ClustalO and extracting the positions of the conserved and unconserved regions
f=open('Allignment.txt')
lines=f.readlines()
Stars=""
myString = lines[17]
myString.index('*')
lines=lines[3:]
# Stars array contains all the stars in the allignment file 
for i in range(14,len(lines),16):
    Stars+=lines[i][22:]


# In[ ]:

# funcion that takes the sequence name and the dictionary containing the data and returns the sequence itself 
def getSequence(a,Dictionary):
    for i in range(len(Dictionary["Sequence Name"])):
        if Dictionary["Sequence Name"][i]==a:
            return Dictionary["Sequences"][i]


# In[ ]:

#getting the refrence sequnce from the dictionary
Refrence=getSequence("NC_045512.2(REF)",Dictionary)




# In[ ]:

#Processing on the conserved Stars array to compute the position of gaps(unconserved) between the conserved regions
import re
newStars=Stars.replace("\n","") 
spaces = re.finditer(r" ", newStars)
idxspaces=[]
for match in spaces:
    run_start = match.start()
    run_end = match.end()
    idxspaces.append(run_end)




# In[ ]:

#computing the unconserved regions and setting them into an array of arrays called uncoservedregions .

j=0
numbr=[]
uncoservedregions=[]
for i in range(len(idxspaces)-1):
    numbr.append(idxspaces[i])
    if idxspaces[i+1]!=idxspaces[i]+1:
        i=i+1
        uncoservedregions.append(numbr)
        numbr=[]




# In[ ]:

#computing the conserved regions and setting them into an array of arrays called coservedregions .

spaces = re.finditer(r"[*]", newStars)
idxstars=[]
for match in spaces:
    run_start = match.start()
    run_end = match.end()
    idxstars.append(run_end)

j2=0
numbr2=[]
coservedregions=[]
for i in range(len(idxstars)-1):
    numbr2.append(idxstars[i])
    if idxstars[i+1]!=idxstars[i]+1:
        i=i+1
        coservedregions.append(numbr2)
        numbr2=[]





# In[ ]:

#reading the reference GFF3 manual 
sheet = pd.read_csv("Refsequence.csv")




# In[ ]:

#extracting the start and end regions and the gene id from the GFF3 manual 
startrange=sheet['Unnamed: 3']
endrange=sheet['Unnamed: 4']
geneId=sheet['Unnamed: 8']



# In[ ]:

# comparing the uncoservedregions with the start and end and matching the regions to each gene 
number3=[]
uncoservedmatchedregions=[]
res = True
count=len(uncoservedregions)
c=1
for ele in uncoservedregions:
    for i in range(2,len(startrange)):
        
        if ele[0] >= startrange[i] and ele[len(ele)-1] <= endrange[i] and ele[0]!=ele[len(ele)-1]:
            if c==1:
                    string="Range="+str(ele[0])+"-"+str(ele[len(ele)-1])
                    number3.append(string)
                    c=0
            number3.append(geneId[i])
            
    if number3:
        uncoservedmatchedregions.append(number3)
    number3=[]
    c=1

            
            




# In[ ]:

#making a dataframe from the array computed
genesunconserved = pd.DataFrame(uncoservedmatchedregions)



# In[ ]:

#saving the genesunconserved to a csv 
genesunconserved.to_csv('Nonconserved Regions.csv')


# In[ ]:

# comparing the conservedregions with the start and end and matching the regions to each gene and setting them into an array 
number4=[]
coservedmatchedregions=[]
res = True
count=len(uncoservedregions)
c=1
for ele in coservedregions:
    for i in range(2,len(startrange)):
        if ele[0] >= startrange[i] and ele[len(ele)-1] <= endrange[i] : 
                if c==1:
                    if ele[0]==ele[len(ele)-1]:
                        number4.append("Range="+str(ele[0]))
                        c=0
                    else:
                        string="Range="+str(ele[0])+"-"+str(ele[len(ele)-1])
                        number4.append(string)
                        c=0
                number4.append(geneId[i])
        else:
             if c==1:
                    if ele[0]==ele[len(ele)-1]:
                        number4.append("Range="+str(ele[0]))
                        c=0
                    else:
                        string="Range="+str(ele[0])+"-"+str(ele[len(ele)-1])
                        number4.append(string)
                        c=0
            
    coservedmatchedregions.append(number4)
    number4=[]
    c=1



# In[ ]:

#making a dataframe from the array computed
genesconserved = pd.DataFrame(coservedmatchedregions)




# In[ ]:

#saving the conserved genes into a csv
genesconserved.to_csv('Conserved Regions.csv')






