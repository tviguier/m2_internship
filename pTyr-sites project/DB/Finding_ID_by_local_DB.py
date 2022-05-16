"""
FINDING_ID_BY_LOCAL
AUTHOR : THOMAS VIGUIER
AIM : 
- FIND ID FROM FaSeq_ID.csv THANK TO SWISSPROT_DB [Homo sapiens]
"""

import pandas as pd
import numpy as np
import re
import time

print("")
print("###########################################################################################")
print("###########################################################################################")
print("")
print("########## I : MERGING DB HOMO SAPIENS & FASEQ FILE")
print("")
print("###########################################################################################")
print("###########################################################################################")
print("")

Homo_sapiens_DB_NAME = "FaSeq_DB_Homo_sapiens_20431.csv"
FaSeq_FILE_NAME = "FaSeq_56001801068863-82113.csv"

df_Homo_sapiens_DB = pd.read_csv(Homo_sapiens_DB_NAME, delimiter = '\t')
df_FaSeq_FILE = pd.read_csv(FaSeq_FILE_NAME, delimiter = '\t')

Stime = time.time()

print("")
print("###########################################################################################")
print(f"##### STEP I.1 : MAIN INFORMATIONS ")
print("###########################################################################################")
print("")

OUTPUT_FORMAT = ".csv"

print(f" - '{Homo_sapiens_DB_NAME}' SELECTED")
print(f" - '{FaSeq_FILE_NAME}' SELECTED")
print(f" - FORMAT : {OUTPUT_FORMAT}")
print("")

print(f" - LENGTH OF {Homo_sapiens_DB_NAME} : {len(df_Homo_sapiens_DB)}")
print(f" - LENGTH OF {FaSeq_FILE_NAME} : {len(df_FaSeq_FILE)}")
print("")

df_Homo_sapiens_DB = df_Homo_sapiens_DB.rename(columns = {'ID': 'HS_ID', 'ID_NAME': 'HS_ID_NAME', 'LEN_SEQ': 'HS_LEN_SEQ', 'SEQUENCE': 'HS_SEQUENCE'})
df_FaSeq_FILE = df_FaSeq_FILE.rename(columns = {'ID': 'FASEQ_ID', 'p':'FASEQ_p', 'LEN_SEQ': 'FASEQ_LEN_SEQ', 'SEQUENCE': 'FASEQ_SEQUENCE'})

print(df_Homo_sapiens_DB)
print("")
print("-----------------------------")
print("")
print(df_FaSeq_FILE)

MERGE_VERSIONS = {0 : ['HS_LEN_SEQ', 'FASEQ_LEN_SEQ'], # MERGE WITH LEN_SEQ
                  1 : ['HS_SEQUENCE', 'FASEQ_SEQUENCE'], # MERGE WITH SEQUENCE
                  2 : ['HS_LEN_SEQ', 'FASEQ_LEN_SEQ', 'HS_SEQUENCE', 'FASEQ_SEQUENCE']} # MERGE WITH ADDING RESPECTIVE LEN_SEQ & SEQUENCE

Homo_sapiens_DB_VERSIONS = {0 : 'HS_LEN_SEQ', # MERGE WITH LEN_SEQ
                            1 : 'HS_SEQUENCE', # MERGE WITH SEQUENCE
                            2 : ['HS_LEN_SEQ', 'HS_SEQUENCE']} # MERGE WITH ADDING RESPECTIVE LEN_SEQ & SEQUENCE


FaSeq_FILE_VERSIONS = {0 : 'FASEQ_LEN_SEQ', # MERGE WITH LEN_SEQ
                       1 : 'FASEQ_SEQUENCE', # MERGE WITH SEQUENCE
                       2 : ['FASEQ_LEN_SEQ', 'FASEQ_SEQUENCE']} # MERGE WITH ADDING RESPECTIVE LEN_SEQ & SEQUENCE


print(f" - {len(MERGE_VERSIONS)} MERGING VERSIONS IDENTIFIED :")
print("")

i = 0 

while i < len(MERGE_VERSIONS) :
    
    print(f"VERSION #{i} : {', '.join(MERGE_VERSIONS[i])}")
    i += 1

print("")
print("")

i = 0

while i < len(MERGE_VERSIONS) :
    
    print(f"##### MERGE VERSION #{i}")
    print("############################################################")
    
    # DATAFRAME CREATION BY MERGING df_REF & df_ALT
    df_MERGE = pd.DataFrame(df_Homo_sapiens_DB.merge(df_FaSeq_FILE, left_on = Homo_sapiens_DB_VERSIONS[i], right_on = FaSeq_FILE_VERSIONS[i] , how = 'outer', indicator = True))
    
    # CHANGING VALUES INSIDE 'PRESENCE' COLUMN
    df_MERGE._merge = ((df_MERGE._merge.str.replace('left_only', 'HS')).replace('right_only', 'FASEQ')).replace('both', 'BOTH')
    
    # COLUMN CREATION '△_PTMscore' TO VISUALIZE THE VALUES DIFFERENCES BETWEEN REF & ALT
    df_MERGE['LEN_SEQ_MATCH'] = np.where(df_MERGE['HS_LEN_SEQ'] == df_MERGE['FASEQ_LEN_SEQ'], 'IDENTICAL', 'DIFFRENT')
    df_MERGE['SEQUENCE_MATCH'] = np.where(df_MERGE['HS_SEQUENCE'] == df_MERGE['FASEQ_SEQUENCE'], 'IDENTICAL', 'DIFFRENT')
    
    # COLUMNS RE-ORDER
    df_MERGE = df_MERGE[['FASEQ_ID', 'FASEQ_p', 'HS_ID', 'HS_ID_NAME', 'FASEQ_TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'HS_LEN_SEQ', 'FASEQ_LEN_SEQ', 'HS_SEQUENCE', 'FASEQ_SEQUENCE', '_merge', 'LEN_SEQ_MATCH', 'SEQUENCE_MATCH']]
    
    # RENAME '_merge' COLUMN BY COLUMN
    MERGE_COLUMN = f"MERGE#{i}_STATUS"
    df_MERGE = df_MERGE.rename(columns={"_merge": MERGE_COLUMN})
    
    #ADD NEW VARIABLES 
    BOTH = df_MERGE[df_MERGE[MERGE_COLUMN] == 'BOTH']
    HS = df_MERGE[df_MERGE[MERGE_COLUMN] == 'HS']
    FASEQ = df_MERGE[df_MERGE[MERGE_COLUMN] == 'FASEQ']
    
    OUTPUT_MERGE_NAME = f"HS-FASEQ_MERGE#{i}"
    OUTPUT_MERGE_FILE = f"{OUTPUT_MERGE_NAME}-{len(BOTH)}-{len(HS)}-{len(FASEQ)}{OUTPUT_FORMAT}"
    
    df_MERGE.to_csv(OUTPUT_MERGE_FILE, sep = '\t', mode='a', index=False)
    
    print("")
    print(f"##### MERGING#{i} RESULTS")
    print("")
    print(f" - NUMBER OF SAME ROW (BOTH) : {len(BOTH)}")
    print(f" - NUMBER OF ROWS ONLY CONTAINED IN {Homo_sapiens_DB_NAME} : {len(HS)}")
    print(f" - NUMBER OF ROWS ONLY CONTAINED IN {FaSeq_FILE_NAME} : {len(FASEQ)}")
    print("")
    print(f" - FILE '{OUTPUT_MERGE_FILE}' CREATED")
    print("")
    print("")

    i += 1

print("############################################################")
print("##### DONE")
print("############################################################")
print("")