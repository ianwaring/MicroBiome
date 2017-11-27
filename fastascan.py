#!/usr/bin/python
"""
fastascan() - Coursera Python for Genomic Science Final Project
Ian Waring, Software Enabled Services Ltd, Monday 28th Sep 2015

Usage: fastascan('Filename')

where Filename is the FASTA format file to be analysed
"""

import sys
import operator

# First some support functions

def has_stop_codon(dna,frame=0) :
    """
    This function checks if a given DNA sequence contains in frame stop codons

    dna is the dna sequence to be analysed
    frame is the offset (0,1 or 2) from start of the sequence to scan triplets
    """
    stop_codons-['tga','tag','taa']
    stop_codon_found = False
    for i in range(frame, len(dna),3):
        codon=dna[i:i+3].lower()
        if codon in stop_codons :
            stop_codon_found=True
            break
        return stop_codon_found

def orfs(dna,frame=0) :
    """
    This function outputs a list of all ORFs found in string 'dna', using
    triplet boundaries defined by frame position 'frame'

    dna is the dna sequence to be analysed
    frame is the offset (0,1 or 2) from start of the sequence to scan triplets
    """
    orfs=[]
    orf=''
    start_codons=['atg']
    stop_codons=['tga','tag','taa']
    start_codon_active = False
    for i in range(frame, len(dna),3):
        codon=dna[i:i+3].lower()
        if start_codon_active == False and codon in start_codons :
            start_codon_active = True
        if start_codon_active and codon in stop_codons :
            orfs.append(orf+codon)
            orf=''
            start_codon_active = False
        else:
            if start_codon_active:
                orf=orf+codon
    return orfs

def complememt(dna) :
    """
    This function returns the reverse complement of a DNA sequence string

    dna is the dna sequence to be translated. Assumes already in lower case
    """
    basecomplement={'a':'t','c':'g','g':'c','t':'a','n':'n'}
    letters = list(dna[::-1])   #blow the string backwards into a letters list
    letters = [basecomplement[base] for base in letters]  #convert all
    return ''.join(letters)     #turn letters back to string

    
filename = 'dna3.fasta'
try:
    f=open(filename,'r')
except IOError:
    print ("File %s does not exist!!!" % filename)



seqs={}     #seqs will contain our dictionary of sequences read in
num_headers=0

for line in f:
    line=line.rstrip()  #remove carriage control characters
    # Test for Header
    if line[0] == ">":
        words=line.split()
        name=words[0][1:]
        seqs[name]=''
        num_headers+=1
    else:
        # Sequence, not a header
        seqs[name]=seqs[name]+line.lower()
f.close()

# Question 1:

print("\nQ1: Number of records in file is %d \n" % num_headers)

# Question 2:

print "\nQuestion 2: Lengths of sequences in file, longest, shortest + ids\n"
i2 = 0
for k in sorted(seqs, key=lambda k: len(seqs[k]), reverse=False):
    i2 += 1
    print i2, k, len(seqs[k])

# Question 3: Length of the Largest ORF in the file, and which Identifier its in
# Question 3: For each identifier, longest length ORF in place
# Question 3: For each identifier, position of the largest ORF inside it

print "\nQuestion 3: Longest ORFs per id per frame with ids\n"

i2 = 0
for k in seqs:
    i2+=1
    for f in range(0,3):
        allorfs=orfs(seqs[k],f)
        if len(allorfs)>0:
            longestorf=max(allorfs, key=len)
            longestorflen=len(longestorf)
            longestorfpos=seqs[k].find(longestorf)+1
            print i2, f+1, k, longestorflen, longestorfpos

# Question 4: find all repeats of supplied length in FASTA file sequences

print "\nQuestion 4: all repeats in Fasta file of specified length\n"

length = int(raw_input("Length of Repeats to find: "))

repeats={}
for k in seqs:
    seq=seqs[k]
    for i in range(0,len(seq)-1):
        item=seq[i:i+length]
        if len(item)<>length:   #Ignore any sub-size superfluos content at end
            pass
        else:
            if item in repeats:
                repeats[item]+=1
            else:
                repeats[item]=1

mostrepeatedkey=max(repeats.iteritems(), key=operator.itemgetter(1))[0]
print "Most Repeated Sequence is:",mostrepeatedkey, "at ",repeats[mostrepeatedkey],"times\n"
for k in repeats:
    if repeats[k] == 8:
        print k,repeats[k]

print "End of Project :-)\n"


