import pandas as pd
import numpy as np
import re
import time

pTyr_FILE = input(">>> INPUT pTyr FILE : ") 
#pTyr_FILE = "pTyr_DB_56001801068863.txt"

print("")

#COMMON_ID = input(f">>> INPUT ID OF '{pTyr_FILE}' : ")
COMMON_ID = "56001801068863"

print("")

DataBase = input(f">>> INPUT BLAST PROTEIN DATABASE (.csv) : ")
#DataBase = "Blast_DB_6038.csv"

POS_INPUTS = ['YES', 'Y', 'yes', 'y', 'oui', 'OUI']

print("")
print(f"##### FILE '{pTyr_FILE}' SELECTED")
print(f"##### FILE '{DataBase}' SELECTED")
print("")

df_pTyr_FILE = pd.read_csv(pTyr_FILE, delimiter = '\t')
pTyr_FORMAT = f".{pTyr_FILE.split('.')[-1]}"
pTyr_FILE_NAME = pTyr_FILE.replace(f'{pTyr_FORMAT}', '')

print(" - pTyr FILE NAME : ", pTyr_FILE_NAME)
print(f" - FORMAT : {pTyr_FORMAT}")
print(f" - NUMBER OF pTyr-Sites DETECTED : {len(df_pTyr_FILE)}")
print("")

df_DataBase = pd.read_csv(DataBase, delimiter = '\t')
Blast_FORMAT = f".{DataBase.split('.')[-1]}"
Blast_FILE_NAME = DataBase.replace(f'{Blast_FORMAT}', '')

print(" - Blast DATABASE FILE NAME : ", Blast_FILE_NAME)
print(f" - FORMAT : {Blast_FORMAT}")
print(f" - NUMBER OF IDs DETECTED : {len(df_DataBase)}")
print("")

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################

print("")
print("RESTRUCTURATION START... (TXT TO DATAFRAME)")
print("")

Stime = time.time()

df_pTyr_FILE = pd.read_csv(pTyr_FILE, delimiter = '\t')

df_pTyr_FILE.columns =['ID', 'POSITION', 'RESIDUE', 'PTMscore']

df_pTyr_RESTRUC = pd.DataFrame()

df_pTyr_RESTRUC['ID'] = np.where((df_pTyr_FILE['ID'].str.split('.', expand=True)[0]) == '>REF', COMMON_ID , COMMON_ID)
df_pTyr_RESTRUC['pTyr_TYPE'] = np.where((df_pTyr_FILE['ID'].str.split('.', expand=True)[0]) == '>REF', 'REF', 'ALT')
df_pTyr_RESTRUC['p'] = (df_pTyr_FILE['ID'].str.split(' ', expand=True)[0]).str.split('.p', expand=True)[1]
df_pTyr_RESTRUC['TYPE'] = (df_pTyr_FILE['ID'].str.split(' ', expand=True)[1]).str.split(':', expand=True)[1]
df_pTyr_RESTRUC['LENGTH'] = (df_pTyr_FILE['ID'].str.split(' ', expand=True)[2]).str.split(':', expand=True)[1]
df_pTyr_RESTRUC['GC'] = (df_pTyr_FILE['ID'].str.split(' ', expand=True)[3]).str.split(':', expand=True)[1]
df_pTyr_RESTRUC['RANGE'] = ((df_pTyr_FILE['ID'].str.split(' ', expand=True)[4]).str.split(':', expand=True)[1]).str.split('(', expand=True)[0]
df_pTyr_RESTRUC['SIGN'] = ((df_pTyr_FILE['ID'].str.split(' ', expand=True)[4]).str.split('(', expand=True)[1]).str.split(')', expand=True)[0]
df_pTyr_RESTRUC['POSITION'] = df_pTyr_FILE['POSITION']
df_pTyr_RESTRUC['RESIDUE'] = df_pTyr_FILE['RESIDUE']
df_pTyr_RESTRUC['PTMscore'] = df_pTyr_FILE['PTMscore']

OUTPUT_pTyr_FORMAT = ".csv" 

OUTPUT_pTyr_NAME = f"pTyr_{COMMON_ID}-{len(df_pTyr_FILE)}"
OUTPUT_pTyr_FILE = f"{OUTPUT_pTyr_NAME}{OUTPUT_pTyr_FORMAT}"

df_pTyr_RESTRUC.to_csv(OUTPUT_pTyr_FILE, sep = '\t', mode='a', index=False)

Etime = time.time()
Rtime = Etime - Stime

print("RESTRUCTURATION DONE")
print("")
print(f" - RUN TIME : {round(Rtime, 3)}")
print("")
print(f" - FILE '{OUTPUT_pTyr_FILE}' CREATED")
print("")
#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################

print("")
print("COLUMNS ADDS START ...")

Stime = time.time()

df_pTyr_RESTRUC_FILE = pd.read_csv(OUTPUT_pTyr_FILE, delimiter = '\t')

df_pTyr_MERGE = pd.DataFrame(df_pTyr_RESTRUC_FILE.merge(df_DataBase, left_on = ['pTyr_TYPE', 'p'], right_on = ['SEQ_TYPE', 'p'] , how = 'outer', indicator = True))

# CHANGING VALUES INSIDE 'PRESENCE' COLUMN
df_pTyr_MERGE._merge = ((df_pTyr_MERGE._merge.str.replace('left_only', 'pTyr_FILE')).replace('right_only', 'BLAST_db')).replace('both', 'BOTH')

df_pTyr_FULL = pd.DataFrame()
df_pTyr_FINAL = pd.DataFrame()

# COLUMNS RE-ORDER
df_pTyr_FULL = df_pTyr_MERGE[['ID', 'pTyr_TYPE', 'p', 'BLAST_ACCESSION', 'BLAST_ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'PTMscore', '_merge']]

# RENAME '_merge' COLUMN BY COLUMN
MERGE_COLUMN = f"PRESENCE"
df_pTyr_FULL = df_pTyr_FULL.rename(columns={"_merge": MERGE_COLUMN})

#ADD NEW VARIABLES 
pTyr_BOTH = df_pTyr_FULL[df_pTyr_FULL[MERGE_COLUMN] == 'BOTH']
pTyr_pTyr_FILE = df_pTyr_FULL[df_pTyr_FULL[MERGE_COLUMN] == 'pTyr_FILE']
pTyr_Blast_db = df_pTyr_FULL[df_pTyr_FULL[MERGE_COLUMN] == 'BLAST_db']

df_pTyr_FINAL = df_pTyr_FULL[df_pTyr_FULL[MERGE_COLUMN] != 'BLAST_db']
df_pTyr_FINAL = df_pTyr_FINAL[['ID', 'pTyr_TYPE', 'p', 'BLAST_ACCESSION', 'BLAST_ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'PTMscore']]

OUTPUT_pTyr_FINAL_NAME = f"{pTyr_FILE_NAME}_FINAL"
OUTPUT_pTyr_FULL_NAME = f"{pTyr_FILE_NAME}_FULL-{len(pTyr_BOTH)}-{len(pTyr_pTyr_FILE)}-{len(pTyr_Blast_db)}"

OUTPUT_pTyr_FINAL_FILE = f"{OUTPUT_pTyr_FINAL_NAME}{OUTPUT_pTyr_FORMAT}"
OUTPUT_pTyr_FULL_FILE = f"{OUTPUT_pTyr_FULL_NAME}{OUTPUT_pTyr_FORMAT}"

df_pTyr_FINAL.to_csv(OUTPUT_pTyr_FINAL_FILE, sep = '\t', mode='a', index=False)
df_pTyr_FULL.to_csv(OUTPUT_pTyr_FULL_FILE, sep = '\t', mode='a', index=False)

Etime = time.time()
Rtime = Etime - Stime

print("COLUMNS ADDS DONE")
print("")
print(f" - RUN TIME : {round(Rtime, 3)}")
print("")
print(f" - NUMBER OF ROWS WITH ACCESSION IDs & NAMES KNOWN : {len(pTyr_BOTH)}")
print(f" - NUMBER OF ROWS WITH ONLY pTyr IDs : {len(pTyr_pTyr_FILE)}")
print(f" - NUMBER OF BLAST IDs NOT IN '{OUTPUT_pTyr_FILE}' : {len(pTyr_Blast_db)}")
print("")
print(f" - FILE '{OUTPUT_pTyr_FINAL_FILE}' CREATED")
print(f" - FILE '{OUTPUT_pTyr_FULL_FILE}' CREATED")
print("")
print("")
print("END")
