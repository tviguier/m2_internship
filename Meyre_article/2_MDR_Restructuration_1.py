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

I = input(">>> INPUT THE MDR FILE TO CONVERT : ") 

Stime = time.time()

print(f"FILE '{I}' SELECTED")
print("")
print("RUNNING ...")
print("")

###################################################
### PROGRAM
###################################################

I_df = pd.read_csv(I, delimiter = '\t')

### SIZE OF THE INPUT FILE 
I_Size = os.path.getsize(I)

ID_Headers = []

len_ini = len(I_df)

i = 0
nb_Headers = 0

print(f"--- MAIN INFORMATIONS ABOUT {I} ---")
print("")
print(f"- SIZE : {I_Size} Bytes ({round((I_Size/1000000),1)} MB)")
print(f"- {len_ini} ROWS DETECTED")
print("")

df_Final = pd.DataFrame()
df_Final[["ID", "Position", "Residue", "PTMscores", "Cutoff=0.5"]] = ""


while i < len_ini :
    
    ID = I_df.loc[i, 'ID'].replace(">", "")
    PTM = str(I_df['PTMscores'].iloc[i])
    
    if PTM[0] == 'P' :
    	
    	PTM = PTM.replace('Phosphotyrosine:', '')
    	PTM = float(PTM)
    	
    	df_Final = df_Final.append({"ID": I_df["ID"].iloc[i], "Position": I_df["Position"].iloc[i], "Residue": I_df["Residue"].iloc[i], "PTMscores": PTM , "Cutoff=0.5": I_df["Cutoff=0.5"].iloc[i]}, ignore_index = True)
    	   	
    	print(I_df.loc[i])
    	print("")
    
    if (I_df.loc[i, 'ID'])[0] == '>' :
        
        ID_Headers.append(ID)
        I_df = I_df.drop([i])
        
        nb_Headers += 1
      	
    	
    i += 1

I_df.reset_index(inplace=True, drop=True)

print(df_Final)
print("")

"""
while i < len_ini :
    
    ID = I_df.loc[i, 'ID'].replace(">", "")
    
    if (I_df.loc[i, 'ID'])[0] == '>' :
        
        ID_Headers.append(ID)
        I_df = I_df.drop([i])
        
        nb_Headers += 1
    
    else : 
    	
    	PTM = str(I_df['PTMscores'].iloc[i])
    	PTM = PTM.replace('Phosphotyrosine:', '')
    	PTM = float(PTM)
        
    	#I_df["PTMscores"].replace({I_df["PTMscores"].loc[i] : PTM}, inplace = True)
    	I_df.replace({I_df["PTMscores"].loc[i] : PTM}, inplace = True)
    	
    	print(I_df.loc[i])
    	print("")
    	
    i += 1

I_df.reset_index(inplace=True, drop=True)

print(I_df)
print("")

i = 0

while i < len(I_df) :
    
    PTM = str(I_df['PTMscores'].iloc[i])
    PTM = PTM.replace('Phosphotyrosine:', '')
    PTM = float(PTM)
    
    print(i, " ", PTM)
    
    #I_df["PTMscores"].replace({I_df["PTMscores"].loc[i] : PTM}, inplace = True)
    I_df.replace({I_df["PTMscores"].loc[i] : PTM}, inplace = True)
    
    i += 1
   
print(I_df)
print("")
print("")

### UNIQUE PROTEIN ID 
I_ID_Unique = list(pd.unique(I_df["ID"].values.ravel()))
ID_Headers_Unique = list(set(ID_Headers))

Etime = time.time()
Rtime = Etime - Stime

print(f"--- MAIN INFORMATIONS ABOUT {I} ---")
print("")
print(f"- SIZE : {I_Size} Bytes ({round((I_Size/1000000),1)} MB)")
print(f"- {len_ini} ROWS DETECTED")
print("")
print("-- SEARCH OF TYROSINE PHOSPHORYLATIONS --")
print("")
print(f"-- {nb_Headers} FASTA HEADERS DETECTED ({len(ID_Headers_Unique)} UNIQUE FASTA HEADERS)")
print(f"-- {len(I_df)} TYROSINE PHOSPHORYLATIONS DETECTED FOR {len(I_ID_Unique)} UNIQUE PROTEIN ID")
print(f"--- {len(ID_Headers_Unique) - len(I_ID_Unique)} FASTA SEQUENCES DOESN'T CONTAIN TYROSINE PHOSPHORYLATIONS")
print("")
print("RUN TIME : ", round(Rtime, 3), "s")
print("----------------------------------------")
print("")

if len(I_ID_Unique) != len(ID_Headers_Unique) :

    ### CREATION OF A TEMPORARY DF WITH MSSING PROTEIN ID (MID)
    df_M_ID = pd.DataFrame()
    df_M_ID["ID"] = ""
    
    i = 0
    
    while i < len(ID_Headers_Unique) :
        
        ID = ID_Headers_Unique[i]
               
        if ID not in I_ID_Unique :
            
            df_M_ID = df_M_ID.append({"ID": ID}, ignore_index = True)
            
        i += 1
    
    print("")
    Rep = input(">>> PRINT THE PROTEINS ID WITH NO TYROSINE PHOSHORYLATIONS ? (Y/N) :")
    print("")
    
    if Rep in ['Y', 'y', 'YES', 'yes']:
           
        print("- DATAFRAME OF PROTEIN ID WITH NO TYROSINE PHOSPHORYLATIONS")
        print("")
        print(df_M_ID)
        print("")

        Rep = str(input(">>> SAVE THE DATAFRAME ? (Y/N) :"))
        print("")

        if Rep in ['Y', 'y', 'YES', 'yes']:

            O_File = "Missing_ID_from_" + I.replace('.txt', '.csv')
            df_M_ID.to_csv(O_File, sep="\t", mode='a', index=True)
            print(f"FILE '{O_File}' CREATED")
            print("")

O_File = I.replace(".txt", ".csv")
I_df.to_csv(O_File, sep="\t", mode='a', index=False)

print(f"FILE '{O_File}' CREATED")
"""
print("")
print("")
print("##### END #####")
