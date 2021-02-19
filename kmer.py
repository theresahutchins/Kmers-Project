# reference video https://www.youtube.com/watch?v=dQG4-Gwo4BE&ab_channel=NextGenGenomics
# file names are 3UTR.fasta and 5UTR.fasta

#packages 
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

		if kmer in kFreq:
			kFreq[kmer] += 1 #if it matches add 1 to frequency 
		else:
			kFreq[kmer] = 1 #if it doesnt match, adds key to dict w/ value of 1

	return kFreq

def runPerKmer():
    finalDic5={}
    finalDic3={}
    kval = int(input("What length of k: "))
    UTR = input("5' or 3' UTR (enter 5 or 3): ")
    if UTR == "5":
        for key in UTR5_Dict:
            finalDic5[key]= kmers(UTR5_Dict[key],kval)
        print(finalDic5)
    elif UTR == "3":
        for key in UTR3_Dict:
            finalDic3[key]= kmers(UTR3_Dict[key],kval)
        print(finalDic3)
    
	#ListOfKmers = [[kmersLengthk[sequence], sequence] for sequence in kmersLengthk] #list the diff kmers
	#ListOfKmers.sort() #sort list so most frequent show up first
	#print(ListOfKmers)
runPerKmer()

#average the five numbers for each time given per gene in white2017_tpm.txt and multiply by frequencies we produced of kmers
#