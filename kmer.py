# reference video https://www.youtube.com/watch?v=dQG4-Gwo4BE&ab_channel=NextGenGenomics
# file names are 3UTR.fasta and 5UTR.fasta

#packages 
import pandas as pd
import numpy as np
from Bio import SeqIO

#open/read/format 5' and 3' fa files
UTR5_Dict={}
with open("5UTR.fasta", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        UTR5_Dict[record.id]= record.seq

UTR3_Dict={}
with open("3UTR.fasta", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        UTR3_Dict[record.id]= record.seq
        
K=0
UTR=0        
#data = pd.read_csv('white2017_tpm.txt', sep="\t", index_col='GeneID')
#data2= pd.DataFrame(index = data.index)
#targetFile1 = open("3UTR.fasta") 
#targetFile2 = open("5UTR.fasta")
#readSeq1 = targetFile1.read()
#readSeq2 = targetFile2.read()

# I Dont think these are formatted corerctly!! 
# not sure how atm (results print some non genetic info)
#FivePrime= "".join(readSeq1.split())
#ThreePrime= "".join(readSeq2.split())

#finding all possible kmers of len k
def kmers(Tseq, k):
    kFreq = {}
    for i in range(0, len(Tseq) - k +1): #range 0 -> length of target seq 
        kmer = Tseq[i:i + k] #indexing way to slice up target sequence according to k
        if "N" not in kmer:
            if kmer in kFreq:
                kFreq[kmer] += 1 #if it matches add 1 to frequency 
            else:
                kFreq[kmer] = 1 #if it doesnt match, adds key to dict w/ value of 1
    return kFreq

def runPerKmer():
    finalDic5={}
    finalDic3={}
    kval = int(input("What length of k: "))
    global K
    K=kval
    utr = input("5' or 3' UTR (enter 5 or 3): ")
    global UTR
    UTR=utr
    if UTR == "5":
        for key in UTR5_Dict:
            finalDic5[key]= kmers(UTR5_Dict[key],kval)
        #print(finalDic5)
        return finalDic5
    elif UTR == "3":
        for key in UTR3_Dict:
            finalDic3[key]= kmers(UTR3_Dict[key],kval)
        #print(finalDic3)
        return finalDic3
	#ListOfKmers = [[kmersLengthk[sequence], sequence] for sequence in kmersLengthk] #list the diff kmers
	#ListOfKmers.sort() #sort list so most frequent show up first
	#print(ListOfKmers)
#runPerKmer()

#average the five numbers for each time given per gene in white2017_tpm.txt, take log base 2 of that, add 0.5, and multiply by frequencies we produced of kmers
def rnaSeq():
    data_white2017 = pd.read_csv('white2017_tpm.txt', sep="\t", index_col='GeneID')
    data_log2= pd.DataFrame(index = data_white2017.index)
    
    ## Get averages for each devolopmental period. Add 0.5 and take log base 2 of that. Put into new dataframe. 
    i=0
    size = len(data_white2017.columns)
    devolpmental_times=[]
    while i < size:
        devolpmental_times.append(str(data_white2017.columns[i])[0:-2])
        name= 'log2(mean of '+ str(data_white2017.columns[i])[0:-2] + ' cells + 0.5)'
        data_log2[name] = np.log2(data_white2017.iloc[:, i:i+5].mean(axis=1)+0.5)
        i=i+5
    ## Rename all columns of data_log2 to match that of data_white2017
    temp = {}
    for i in data_log2.index:
        r = 'danRer11_ensGene_' + i
        temp[i]=r
    data_log2=data_log2.rename(temp)
    
    ## Get kmer counts by runing runPerKmer and put into dataframe
    freqs= runPerKmer()
    df_freqs = pd.DataFrame.from_dict(freqs, orient='index')
    
    #run for specific time 
    time_name='log2(mean of 1-cell cells + 0.5)' ### change the name from the specific time you want 
    time1= data_log2[time_name] 
    # change all NA to 0
    
    df_freqs.fillna(value=0,inplace=True) 
    finDict={}
    #time= data2.columns
    for i in data_log2.index:
        if i in df_freqs.index:
            df_freqs.loc[i] = df_freqs.loc[i]*time1[i]
    # Get rid of kmers containining N
    cols = [c for c in df_freqs.columns if 'N' not in c ]
    df_freqs= df_freqs[cols]
    
    output_name= str(UTR) + "'UTR_" +str(K)+"mer_"+ devolpmental_times[0]+ ".txt" #will need to change devolpmental_times[0] when you run for the other times.
    df_freqs.to_csv(output_name, sep= "\t")    
rnaSeq()