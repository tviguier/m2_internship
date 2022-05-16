
# SPLITING PART

##################################################
### IMPORT PACKAGES
###################################################

import pandas as pd
from Bio import SeqIO
import os
import time

###################################################
### INPUTS
###################################################

I = input(">>> INPUT THE FASTA FILE TO SPLIT : ")

Stime = time.time()

print(f"FILE '{I}' SELECTED")
print("")
print("RUNNING ...")

###################################################
### PROGRAM
###################################################

I_Size = os.path.getsize(I) #SIZE OF THE INPUT FILE 

I_FaSeqs = list(SeqIO.parse(I, 'fasta'))

Etime = time.time()
Rtime = Etime - Stime

### PRINT TO USER
print("")
print(f"##### MAIN INFORMATIONS of {I} #####")
print("")
print(f"- SIZE : {I_Size} Bytes ({round((I_Size/1000000),1)} MB)")
print(f"- {len(I_FaSeqs)} FASTA SEQUENCES DETECTED")
print("")

I_part = int(input(">>> NUMBER OF PART WANTED :"))
print("")

I_rows = int(round(len(I_FaSeqs)/I_part, 0))

print("PART 1 TO {} : {} FASTA SEQUENCES".format(I_part - 1, I_rows))
print("PART {} : {} FASTA SEQUENCES".format(I_part, len(I_FaSeqs) - (I_rows * (I_part - 1))))
print("")

i = 0
j = 0

seq_start = 0
seq_stop = 0

while i < I_part - 1:
    
    O_File_fa = I.replace(".fa", "_part{}.fa".format(i+1))
    
    O_File = open(O_File_fa, "w+")

    seq_stop = seq_start + I_rows

    print(f"PART N°{i + 1} WITH {seq_stop - seq_start} FASTA SEQUENCES")
    print(f"- FILE '{O_File_fa}' CREATED")
    print("")
    
    for j in range(seq_start, seq_stop):
               
        O_File.write(">" + I_FaSeqs[j].id + "\n")
        O_File.write("".join(list(I_FaSeqs[j].seq)) + "\n")
        O_File.write("\n")
        
    seq_start = seq_stop
    
    j = seq_start
    
    i += 1
    
    O_File.close()

O_File_fa = I.replace(".fa", "_part{}.fa".format(I_part))

O_File = open(O_File_fa, "w+")


print(f"PART N°{ I_part} WITH {len(I_FaSeqs) - seq_start} FASTA SEQUENCES")
print(f"- FILE '{O_File_fa}' CREATED")
print("")

for j in range(seq_start, len(I_FaSeqs)):
                   
        O_File.write(">" + I_FaSeqs[j].id + "\n")
        O_File.write("".join(list(I_FaSeqs[j].seq)) + "\n")
        O_File.write("\n")

seq_start = seq_stop

j = seq_start

i += 1

O_File.close()

print("##### END #####")
print("")
print("!!! WARNING !!!")
print("ALL THESE FILES MUST BE ANALYZE BY MUSITEDEEP PREDICTION ALGO BEFORE TO RUN THE SECOND PYTHON SCRIPT")
print("")