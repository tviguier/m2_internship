
# INFORMATIONS PART

###################################################
### IMPORT PACKAGES
###################################################

import pandas as pd
from Bio import SeqIO
import os
import time

###################################################
### INPUTS
###################################################

I = input(">>> INPUT THE FASTA ANALYZE : ")

Stime = time.time()

print(f"FILE '{I}' SELECTED")
print("")
print("RUNNING ...")

###################################################
### PROGRAM
###################################################

I_Size = os.path.getsize(I) #SIZE OF THE INPUT FILE 

I_FaSeqs = list(SeqIO.parse(I, 'fasta'))

### DATAFRAME CREATION TO STOCK INFORMATIONS FOR EACH FASTA SEQUENCE
I_df = pd.DataFrame()
I_df[["Uniprot_ID", "Name", "Description", "Sequence"]] = ""

i = 0

O_File_fa = str(len(I_FaSeqs)) + "_FaSeqs.fa"
O_File_csv = str(len(I_FaSeqs)) + "_FaSeqs.csv"

O_File = open(O_File_fa, "w+")

while i < len(I_FaSeqs) :
    
    ID = str((I_FaSeqs[i].id).split('|')[1])
    NAME = str((I_FaSeqs[i].id).split('|')[2])
    DESCRIPTION = str(I_FaSeqs[i].description)
    SEQUENCE = "".join(list(I_FaSeqs[i].seq))
    
    I_df = I_df.append({"Uniprot_ID": ID, "Name": NAME, "Description": DESCRIPTION, "Sequence": SEQUENCE}, ignore_index = True)
    
    O_File.write(">" + ID + "\n")
    O_File.write(SEQUENCE + "\n")
    O_File.write("\n")
    
    i += 1

### CREATION OF .csv & .fa FILES
O_File.close()
I_df.to_csv(O_File_csv, sep="\t", mode='a', index=False)


### LIST OF UNIQUE INFORMATIONS INSIDE DATAFRAME --> U_xxx
U_ID = list(pd.unique(I_df["Uniprot_ID"].values.ravel()))
U_NAME = list(pd.unique(I_df["Name"].values.ravel()))
U_DESCRIPTION = list(pd.unique(I_df["Description"].values.ravel()))
U_SEQUENCE = list(pd.unique(I_df["Sequence"].values.ravel()))

Etime = time.time()
Rtime = Etime - Stime

### PRINT TO USER
print("")
print(f"##### MAIN INFORMATIONS of {I} #####")
print("")
print(f"- SIZE : {I_Size} Bytes ({round((I_Size/1000000),1)} MB)")
print(f"- {len(I_FaSeqs)} FASTA SEQUENCES DETECTED")
print("")
print("- FORMATATION -")
print("")
print(f"- 2 FILES HAVE BEEN CREATED :")
print(f"-- '{O_File_fa}'")
print(f"-- '{O_File_csv}'")
print("")
print("-- SEARCH OF DUPLICATES --")
print("")
print(f"-- {len(U_ID)} UNIQUE FASTA SEQUENCES ACCORDING THEIR ID ({len(I_FaSeqs) - len(U_ID)} DUPLICATES)")
print(f"-- {len(U_NAME)} UNIQUE FASTA SEQUENCES ACCORDING THEIR NAME ({len(I_FaSeqs) - len(U_NAME)} DUPLICATES)")
print(f"-- {len(U_DESCRIPTION)} UNIQUE FASTA SEQUENCES ACCORDING THEIR DESCRIPTION ({len(I_FaSeqs) - len(U_DESCRIPTION)} DUPLICATES)")
print(f"-- {len(U_SEQUENCE)} UNIQUE FASTA SEQUENCES ACCORDING THEIR SEQUENCE ({len(I_FaSeqs) - len(U_SEQUENCE)} DUPLICATES)")
print("")
print("RUN TIME : ", round(Rtime, 3), "s")
print("----------------------------------------")
print("")

if len(U_SEQUENCE) != len(I_FaSeqs) :
    
    ### CREATION OF A TEMPORARY DF WITH DUPLICATE VALUE ACCORDING A SPECIFIC COLUMN
    DRows_df = I_df[I_df.duplicated('Sequence')]
    D_df = pd.DataFrame() # D_df = Duplicate df with all duplicates values
    
    i = 0

    while i < len(DRows_df["Sequence"]) :

        Dtemp_df = I_df.loc[I_df['Sequence'] == DRows_df["Sequence"].iloc[i]]
        D_df = D_df.append(Dtemp_df)

        i += 1
    
    Rep = str(input(">>> PRINT THE DUPLICATES VALUES ? (Y/N) :"))
    print("")
    
    if Rep in ['Y', 'y', 'YES', 'yes']:
           
        print("--- DATAFRAME WITH ALL DUPLICATES ---")
        print("")
        print(D_df[['Sequence', 'Uniprot_ID', 'Name']])
        print("")

        Rep = str(input(">>> SAVE THE DUPLICATE DATAFRAME ? (Y/N) :"))
        print("")

        if Rep in ['Y', 'y', 'YES', 'yes']:

            O_File_D = "Duplicates_from_" + O_File_csv
            D_df.to_csv(O_File_D, sep="\t", mode='a', index=False)
            #D_df[['Sequence', 'Uniprot_ID', 'Name']].to_csv(O_File, sep="\t", mode='a', index=False)
            print(f"FILE '{O_File_D}' CREATED")
            print("")

print("")
print("##### END #####")
