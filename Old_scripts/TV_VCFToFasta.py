from Bio import SeqIO
from art import *
import sys
import io
import pandas as pd
import datetime
import time
import art

tprint("VCF   TO   FASTA")
print("AUTHOR : THOMAS VIGUIER")
print("")

IFile = input("> VCF FILE AS INPUT : ") # file to parse et compare with reference file
RFile = input(">> FASTA FILE AS REFERENCE : ") # genuine file to make comparisons with the mod_file
OFile = input(">>> FASTA FILE AS OUTPUT : ") # "Fa2VCF_out07.vcf" file with written mismatches between mod_file et reference_file

print("")
print("-------------------------------------------")

CurrentDate=datetime.date.today()

Stime = time.time()

RSeq = list(SeqIO.parse(RFile, 'fasta'))

def ReadVCF(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str}, sep='\t').rename(columns={'#CHROM': 'CHROM'})

df = ReadVCF(IFile)

i = 0
RSeqIDs = []

while i < len(RSeq) :
    RSeqIDs.append(RSeq[i].id)
    i += 1

ISeqIDs = df["CHROM"].values.ravel()
ISeqIDs = pd.unique(ISeqIDs)
ISeqIDs = list(ISeqIDs)

# --------------------- DF WITH SNP FROM INPUT VCF --------------------

df1 = pd.DataFrame()

df1['ID'] = ""
df1['Description'] = ""
df1['Seq'] = ""

df1 = df1[['ID', 'Description', 'Seq']]

# ---------------------------------------------------------------------

i = 0
rs = 0 # recognized sequences from input file

while i < len(ISeqIDs) :
    
    if ISeqIDs[i] in RSeqIDs:

        RSeq_index = RSeqIDs.index(ISeqIDs[i]) # find the seq ref position inside the ref file
        df1 = df1.append({'ID': RSeq[RSeq_index].id, 'Description': RSeq[RSeq_index].description, 'Seq': str(RSeq[RSeq_index].seq.upper())}, ignore_index=True)
        rs += 1
        
    else:

        df1 = df1.append({'ID': ISeqIDs[i], 'Description': RSeq[RSeq_index].description, 'Seq': "UNKNOW SEQUENCE"}, ignore_index=True)
    
    i += 1

i = 0
ms = 0 # number of mismatchs

while i < len(df) :
    
    if df.CHROM[i] in ISeqIDs:
        
        ISeq_index = ISeqIDs.index(df.CHROM[i])

        Seq = df1.Seq[ISeq_index]

        pos = df.POS[i]
        val = df.ALT[i]
        tmp = list(Seq)
        tmp[pos] = val
        Seq = "".join(tmp)
    
        df1.iat[int(ISeq_index), 2] = Seq
        
        ms += 1
        
    i += 1

i = 0

OFile = open(OFile, "w")

while i < len(df1) :
        
    OFile.write(">" + df1.Description[i] + "\n")
    OFile.write(df1.Seq[i] + "\n")
    OFile.write("" + "\n")
    
    i += 1

OFile.close()

Etime = time.time()
Rtime = Etime - Stime

print("FASTA REFERENCE FILE :", RFile)
print("-- NUMBER OF SEQUENCES :", len(RSeq))
print("")
print("VCF INPUT FILE :", IFile)
print("-- NUMBER OF SEQUENCES :", len(df1))
print("-- NUMBER OF RECOGNIZED SEQUENCES :", rs)
print("-- TOTAL NUMBER OF MISMATCHS :", ms)
print("-------------------------------------------")
print("FASTA OUTPUT CREATED")
print("RUN TIME :", round(Rtime, 3), "s")
print("")
