from Bio import SeqIO
from art import *
import sys
import pandas as pd
import datetime
import time
import art

tprint("FASTA   TO   VCF")
print("AUTHOR : THOMAS VIGUIER")
print("")

IFile = input("> FASTA FILE AS INPUT : ") # "Input_m10seq.fa" file to parse et compare with reference file
RFile = input(">> FASTA FILE AS REFERENCE : ") # "Test999.fa" genuine file to make comparisons with the mod_file
OFile = input(">>> VCF FILE AS OUTPUT : ") # "Fa2VCF_out07.vcf" file with written mismatches between mod_file et reference_file

print("")
print("-------------------------------------------")

CurrentDate=datetime.date.today()

Stime = time.time()

ISeq = list(SeqIO.parse(IFile, 'fasta'))
RSeq = list(SeqIO.parse(RFile, 'fasta'))

i = 0
RSeqIDs = []

while i < len(RSeq) :
    RSeqIDs.append(RSeq[i].id)
    i += 1

# -------------------- VCF FILE OUTPUT CREATION --------------------

df = pd.DataFrame()

df['#CHROM'] = ""
df['POS'] = ""
df['ID'] = ""
df['REF'] = ""
df['ALT'] = ""
df['QUAL'] = ""
df['FILTER'] = ""
df['INFO'] = ""

df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']] 

# ------------------------------------------------------------------

i = 0
rs = 0 # recognized fasta sequences from input
#tms = 0 # total mismatch detected
ms = 0 # mismatch count 

while i < len(ISeq) :
    
    bp = 0 # nt position

    if ISeq[i].id in RSeqIDs :
        
        RSeq_index = RSeqIDs.index(ISeq[i].id) # find the seq ref position inside the ref file
        
        while bp < len(RSeq[RSeq_index].seq) :
            
            if ISeq[i].seq[bp] != RSeq[RSeq_index].seq[bp]:
                
                df = df.append({'#CHROM': RSeq[RSeq_index].id , 'POS': bp, 'ID': ".", 'REF': RSeq[RSeq_index].seq[bp].upper() , 'ALT': ISeq[i].seq[bp].upper(), 'INFO': RSeq[RSeq_index].description}, ignore_index=True)
                ms += 1
                
            bp += 1
    
        if ISeq[i].id in RSeqIDs and ms == 0 :
        
            df = df.append({'#CHROM': ISeq[i].id , 'POS': 0, 'ID': ".", 'REF': RSeq[RSeq_index].seq[0].upper() , 'ALT': RSeq[RSeq_index].seq[0].upper(), 'INFO': ISeq[i].description}, ignore_index=True)
        
        rs += 1
        #tms += ms
    
    else:
        
        df = df.append({'#CHROM': ISeq[i].id , 'POS': "0", 'ID': "ERROR", 'REF': "ERROR" , 'ALT': "ERROR", 'QUAL': "ERROR", 'FILTER': "ERROR", 'INFO': ISeq[i].description}, ignore_index=True)
        
    i += 1

Etime = time.time()
Rtime = Etime - Stime

header = """##fileformat=VCFv4.1
##fileDate={}
##INFO=<ReferenceFile="{}",Number_of_sequences={}>
##INFO=<InputFile="{}",TotalSequences={},TotalRecognized={},TotalMismatchs={}>
##INFO=<RunTime={}s>
""".format(CurrentDate, RFile, len(RSeq), IFile, len(ISeq), rs, ms, round(Rtime, 3))

with open(OFile, 'w') as vcf:
    vcf.write(header)

df.to_csv(OFile, sep="\t", mode='a', index=False)

print("FASTA REFERENCE FILE :", RFile)
print("-- NUMBER OF SEQUENCES :", len(RSeq))
print("")
print("FASTA INPUT FILE :", IFile)
print("-- NUMBER OF SEQUENCES :", len(ISeq))
print("-- NUMBER OF RECOGNIZED SEQUENCES :", rs)
print("-- TOTAL NUMBER OF MISMATCHS :", ms)
print("-------------------------------------------")     
print(df[['#CHROM', 'POS', 'ID', 'REF', 'ALT']])
print("SIZE :",sys.getsizeof(df), "bytes")
print("-------------------------------------------")
print("VCF OUTPUT CREATED")
print("RUN TIME :", round(Rtime, 3), "s")
print("")