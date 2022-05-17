import pandas as pd
import numpy as np
import time 

#Blast_DataBase = input(">>> INPUT Blast DATABASE : ") 
Blast_DataBase = "partial_ID_56001801068863.txt"

print("")

#COMMON_ID = input(f">>> INPUT ID OF '{Blast_DataBase}' : ")
COMMON_ID = "56001801068863"

print("")

print(f"##### FILE '{Blast_DataBase}' SELECTED")
print("")

df_Blast_DataBase = pd.read_csv(Blast_DataBase, delimiter = '\t')
Blast_FORMAT = f".{Blast_DataBase.split('.')[-1]}"
Blast_DataBase_NAME = Blast_DataBase.replace(f'.{Blast_FORMAT}', '')

print(" - FILE NAME : ", Blast_DataBase_NAME)
print(" - ID : ", COMMON_ID)
print(f" - FORMAT : {Blast_FORMAT}")
print(f" - NUMBER OF Blast IDs DETECTED : {len(df_Blast_DataBase)}")
print("")
print("")
print("DATABASE RESTRUCTURATION ... ")
print("")

Stime = time.time()

df_Blast_temp = df_Blast_DataBase[df_Blast_DataBase.columns[[0,1]]]
df_Blast_temp.columns =['ID', 'BLAST']
df_Blast_temp = df_Blast_temp.drop_duplicates()
df_Blast_temp = df_Blast_temp.reset_index(drop=True)

df_Blast_RESTRUC = pd.DataFrame()

df_Blast_RESTRUC['INDIVIDUAL_ID'] = np.where((df_Blast_temp['ID'].str.split('.', expand=True)[0]) == 'REF', COMMON_ID , COMMON_ID)
df_Blast_RESTRUC['SEQ_TYPE'] = np.where((df_Blast_temp['ID'].str.split('.', expand=True)[0]) == 'REF', 'REF', 'ALT')
df_Blast_RESTRUC['p'] = df_Blast_temp['ID'].str.split('.p', expand=True)[1]
df_Blast_RESTRUC['BLAST_ACCESSION'] = df_Blast_temp['BLAST'].str.split('|', expand=True)[1]
df_Blast_RESTRUC['BLAST_ID'] = (df_Blast_temp['BLAST'].str.split('|', expand=True))[2]

OUTPUT_Blast_FORMAT = ".csv"

OUTPUT_Blast_NAME = f"Blast_db_{COMMON_ID}-{len(df_Blast_RESTRUC)}"
OUTPUT_Blast_DataBase = f"{OUTPUT_Blast_NAME}{OUTPUT_Blast_FORMAT}"

df_Blast_RESTRUC.to_csv(OUTPUT_Blast_DataBase, sep='\t', mode='a', index=False)

Etime = time.time()
Rtime = Etime - Stime

print("Restructured Blast DataBase\n\n", df_Blast_RESTRUC)
print("")
print("DATABASE RESTRUCTURATION DONE")
print("")
print(f" - RUN TIME : {round(Rtime, 3)}")
print("")
print(f" - NUMBER OF UNIQUE IDs : {len(df_Blast_RESTRUC)}")
print(f" - FILE '{OUTPUT_Blast_DataBase}' CREATED")
print("")
print("")
print("DONE")
