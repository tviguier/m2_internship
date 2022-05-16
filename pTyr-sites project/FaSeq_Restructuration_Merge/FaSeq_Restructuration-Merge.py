"""
FaSeq RESTRUCTURATION
AUTHOR : THOMAS VIGUIER
AIM : 
- CONVERT FASTA FILE FILE TO DATAFRAME
- SPLITING FASTA REF & ALT ROWS INTO 2 DATAFRAMES
- MERGE THE DATAFRAMES TO SEE THE COMMON AND DIFFERENT FASTA SEQUENCES
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
import time


POS_INPUTS = ['YES', 'Y', 'yes', 'y', 'oui', 'OUI']

print("")
print("########################################################################################################################")
print("########################################################################################################################")
print("")
print("############### FaSeq ANALYZE")
print("")
print("########################################################################################################################")
print("########################################################################################################################")
print("")

FaSeq_FILE = input(">>> INPUT FASTA FORMAT FILE : ") 

print("")

COMMON_ID = input(f">>> INPUT ID OF '{FaSeq_FILE}' : ")
    
print("")
print("###########################################################################################")
print(f"##### MAIN INFORMATIONS ABOUT {FaSeq_FILE}")
print("###########################################################################################")
print("")

list_FaSeq_FILE = list(SeqIO.parse(FaSeq_FILE, "fasta"))
FaSeq_FORMAT = f".{FaSeq_FILE.split('.')[-1]}"
FaSeq_FILE_NAME = FaSeq_FILE.replace(f'.{FaSeq_FORMAT}', '')

print(" - FILE NAME : ", FaSeq_FILE_NAME)
print(" - ID : ", COMMON_ID)
print(f" - FORMAT : {FaSeq_FORMAT}")
print(f" - NUMBER OF FaSeq-Sites DETECTED : {len(list_FaSeq_FILE)}")
print("")

OUTPUT_FaSeq_FORMAT = ".csv" 

OUTPUT_FaSeq_NAME = f"FaSeq_{COMMON_ID}-{len(list_FaSeq_FILE)}"
OUTPUT_FaSeq_FILE = f"{OUTPUT_FaSeq_NAME}{OUTPUT_FaSeq_FORMAT}"

print("")
print("###########################################################################################")
print(f"##### STEP I.2 : RESTRUCTURATION")
print("###########################################################################################")
print("")

Stime = time.time()

print("RESTRUCTURATION ... (FASTA TO DATAFRAME)")
print("")

df_FaSeq_FILE = pd.DataFrame()

i = 0

while i < 8000 : #len(list_FaSeq_FILE) :
    
    ID = (list_FaSeq_FILE[i].id).split(".")[0]
    FASEQ_TYPE = "REF"
    p = (((list_FaSeq_FILE[i].id).split(".")[1]).split(" ")[0]).replace("p", "")
    TYPE = ((list_FaSeq_FILE[i].description).split(" ")[1]).replace("type:", "")
    LENGTH = ((list_FaSeq_FILE[i].description).split(" ")[2]).replace("len:", "")
    GC = ((list_FaSeq_FILE[i].description).split(" ")[3]).replace("gc:", "")
    RANGE = (((list_FaSeq_FILE[i].description).split(" ")[4]).split(":")[1]).split("(")[0]
    SIGN = (((list_FaSeq_FILE[i].description).split(" ")[4]).split("(")[1]).replace(")", "")
    FASEQ = str(list_FaSeq_FILE[i].seq)
    
    #HEADER = str(f"type:{TYPE} len:{LENGTH} gc:{GC} range:{RANGE}({SIGN})")
    
    if ID == str(COMMON_ID) :
        FASEQ_TYPE = "ALT"

    # df_ROW : ROW WITH ALL THE INFORMATIONS CONTAINS IN THE FASTA SEQUENCE
    df_ROW = pd.DataFrame({"ID": [int(COMMON_ID)],
                           "FASEQ_TYPE": [str(FASEQ_TYPE)],
                           "p": [int(p)], 
                           "TYPE": [str(TYPE)] , 
                           "LENGTH": [int(LENGTH)] , 
                           "GC": [str(GC)], 
                           "RANGE": [str(RANGE)],
                           "SIGN": [str(SIGN)],
                           "LEN_SEQ": [int(len(FASEQ))],
                           "SEQUENCE": [str(FASEQ)]})

    df_FaSeq_FILE = pd.concat([df_FaSeq_FILE, df_ROW])
    
    i += 1

df_FaSeq_FILE = df_FaSeq_FILE.reset_index(drop=True)

Etime = time.time()
Rtime = Etime - Stime

print("RESTRUCTURATION DONE")
print("")
print(f" - RUN TIME : {round(Rtime, 3)}")
print("")

FaSeq_INPUT_0 = input(">>> PRINT THE RESULT DATAFRAME ? (Y/N) : ")

print("")

if FaSeq_INPUT_0 in POS_INPUTS :
    
    print("")
    print(df_FaSeq_FILE)
    print("")

FaSeq_INPUT_1 = input(">>> GROUP THE COLUMNS ['TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN'] IN A NEW DATAFRAME ? (Y/N) : ")

if FaSeq_INPUT_1 in POS_INPUTS :
    
    df_FaSeq_HEADER = pd.DataFrame()
    
    df_FaSeq_HEADER['ID'] = df_FaSeq_FILE['ID']
    df_FaSeq_HEADER['p'] = df_FaSeq_FILE['p']
    df_FaSeq_HEADER['FASEQ_TYPE'] = df_FaSeq_FILE['FASEQ_TYPE']
    df_FaSeq_HEADER['HEADER'] = "type:" + df_FaSeq_FILE['TYPE'] + " len:" + df_FaSeq_FILE['LENGTH'].astype(str) + " gc:" + df_FaSeq_FILE['GC'] + " range:" + df_FaSeq_FILE['RANGE'].astype(str) + "(" + df_FaSeq_FILE['SIGN'] + ")"
    df_FaSeq_HEADER['LEN_SEQ'] = df_FaSeq_FILE['LEN_SEQ']
    df_FaSeq_HEADER['SEQUENCE'] = df_FaSeq_FILE['SEQUENCE']
    
    df_FaSeq_HEADER = df_FaSeq_HEADER[['ID', 'FASEQ_TYPE', 'p', 'HEADER', 'LEN_SEQ', 'SEQUENCE']]
    
    print("")
    print("GROUP DONE")
    print("")
    
    FaSeq_INPUT_2 = input(">>> PRINT THE HEADER DATAFRAME ? (Y/N) : ")
    
    if FaSeq_INPUT_2 in POS_INPUTS :
        
        print("")
        print(df_FaSeq_HEADER)

print("")
print("###########################################################################################")
print(f"##### STEP I.3 : SPLITING FaSeq-REF & FaSeq-ALT ROWS")
print("###########################################################################################")
print("")

print("SPLITING ...")
print("")

Stime = time.time()

df_FaSeq_REF = df_FaSeq_FILE[df_FaSeq_FILE['FASEQ_TYPE'] == 'REF']
df_FaSeq_REF = df_FaSeq_REF.rename(columns={'p' : 'REF_p', 'LEN_SEQ' : 'REF_LEN_SEQ', 'SEQUENCE' : 'REF_SEQUENCE'})
df_FaSeq_REF = df_FaSeq_REF[['ID', 'FASEQ_TYPE', 'REF_p', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_LEN_SEQ', 'REF_SEQUENCE']]
                 
df_FaSeq_ALT = df_FaSeq_FILE[df_FaSeq_FILE['FASEQ_TYPE'] == 'ALT']
df_FaSeq_ALT = df_FaSeq_ALT.rename(columns={'p' : 'ALT_p', 'LEN_SEQ' : 'ALT_LEN_SEQ', 'SEQUENCE' : 'ALT_SEQUENCE'})
df_FaSeq_ALT = df_FaSeq_ALT[['ID', 'FASEQ_TYPE', 'ALT_p', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'ALT_LEN_SEQ', 'ALT_SEQUENCE']]

Etime = time.time()
Rtime = Etime - Stime


print("SPLITING DONE")
print("")
print(f" - RUN TIME : {round(Rtime, 3)}")
print("")

FaSeq_INPUT_3 = input(">>> PRINT FaSeq-REF & FaSeq-ALT DATAFRAMES ? (Y/N) : ")

print("")

if FaSeq_INPUT_3 in POS_INPUTS :
    
    print(" - df_FaSeq_REF\n", df_FaSeq_REF)
    print("")
    print(" - df_FaSeq_ALT\n",df_FaSeq_ALT)

print("")

FaSeq_INPUT_4 = input(">>> SAVE FaSeq_RESULTS DATAFRAMES ? (Y/N) : ")

print("")

OUTPUT_FaSeq_REF_NAME = f"{OUTPUT_FaSeq_NAME}_REF-{len(df_FaSeq_REF)}"
OUTPUT_FaSeq_REF_FILE = f"{OUTPUT_FaSeq_REF_NAME}{OUTPUT_FaSeq_FORMAT}"

OUTPUT_FaSeq_ALT_NAME = f"{OUTPUT_FaSeq_NAME}_ALT-{len(df_FaSeq_ALT)}"
OUTPUT_FaSeq_ALT_FILE = f"{OUTPUT_FaSeq_ALT_NAME}{OUTPUT_FaSeq_FORMAT}"

if FaSeq_INPUT_4 in POS_INPUTS :

    df_FaSeq_FILE.to_csv(OUTPUT_FaSeq_FILE, sep='\t', mode='a', index=False)
    df_FaSeq_REF.to_csv(OUTPUT_FaSeq_REF_FILE, sep='\t', mode='a', index=False)
    df_FaSeq_ALT.to_csv(OUTPUT_FaSeq_ALT_FILE, sep='\t', mode='a', index=False)
    
    print(f" - FILE '{OUTPUT_FaSeq_FILE}' CREATED")
    
    if FaSeq_INPUT_1 in POS_INPUTS :
        
        OUTPUT_FaSeq_HEADER_NAME = f"{OUTPUT_FaSeq_NAME}_HEADER-{len(df_FaSeq_HEADER)}"
        OUTPUT_FaSeq_HEADER_FILE = f"{OUTPUT_FaSeq_HEADER_NAME}{OUTPUT_FaSeq_FORMAT}"
        
        df_FaSeq_HEADER.to_csv(OUTPUT_FaSeq_HEADER_FILE, sep='\t', mode='a', index=False)
        
        print(f" - FILE '{OUTPUT_FaSeq_HEADER_FILE}' CREATED")
          
    print(f" - FILE '{OUTPUT_FaSeq_REF_FILE}' CREATED")
    print(f" - FILE '{OUTPUT_FaSeq_ALT_FILE}' CREATED")
    print("")


print("")
print("###########################################################################################")
print(f"##### FaSeq ANALYZE DATAFRAMES CREATION DONE")
print("###########################################################################################")
print("")
print("")

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


print("")
print("###########################################################################################")
print("###########################################################################################")
print("")
print("########## II : MERGING FaSeq-REF & FaSeq-ALT DATAFRAMES")
print("")
print("###########################################################################################")
print("###########################################################################################")
print("")

df_FaSeq_REF = pd.read_csv(OUTPUT_FaSeq_REF_FILE, delimiter = '\t')
df_FaSeq_ALT = pd.read_csv(OUTPUT_FaSeq_ALT_FILE, delimiter = '\t')

Stime = time.time()

OUTPUT_FaSeq_FORMAT = ".csv"

print(f" - '{OUTPUT_FaSeq_REF_FILE}' SELECTED")
print(f" - '{OUTPUT_FaSeq_ALT_FILE}' SELECTED")
print(f" - FORMAT : {OUTPUT_FaSeq_FORMAT}")
print("")

print(f" - LENGTH OF {OUTPUT_FaSeq_REF_FILE} : {len(df_FaSeq_REF)}")
print(f" - LENGTH OF {OUTPUT_FaSeq_ALT_FILE} : {len(df_FaSeq_ALT)}")

FaSeq_MERGE_VERSIONS = {0 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN'], # COMMON MERGE
                        1 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_p', 'ALT_p'], # MERGE WITH ADDING RESPECTIVE p
                        2 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_LEN_SEQ', 'ALT_LEN_SEQ'], # MERGE WITH ADDING RESPECTIVE PTMscore
                        3 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_SEQUENCE', 'ALT_SEQUENCE'],
                        4 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_p', 'ALT_p', 'REF_LEN_SEQ', 'ALT_LEN_SEQ'],
                        5 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_p', 'ALT_p', 'REF_SEQUENCE', 'ALT_SEQUENCE'],
                        6 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_LEN_SEQ', 'ALT_LEN_SEQ', 'REF_SEQUENCE', 'ALT_SEQUENCE'],
                        7 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_p', 'ALT_p', 'REF_LEN_SEQ', 'ALT_LEN_SEQ','REF_LEN_SEQ', 'ALT_LEN_SEQ', 'REF_SEQUENCE', 'ALT_SEQUENCE']} # MERGE WITH ADDING RESPECTIVE p & PTMscore

FaSeq_REF_VERSIONS = {0 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN'],
                      1 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_p'],
                      2 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_LEN_SEQ'],
                      3 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_SEQUENCE'],
                      4 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_p', 'REF_LEN_SEQ'],
                      5 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_p', 'REF_SEQUENCE'],
                      6 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_LEN_SEQ', 'REF_SEQUENCE'],
                      7 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_p', 'REF_LEN_SEQ', 'REF_SEQUENCE']}

FaSeq_ALT_VERSIONS = {0 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN'],
                      1 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'ALT_p'],
                      2 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'ALT_LEN_SEQ'],
                      3 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'ALT_SEQUENCE'],
                      4 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'ALT_p', 'ALT_LEN_SEQ'],
                      5 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'ALT_p', 'ALT_SEQUENCE'],
                      6 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'ALT_LEN_SEQ', 'ALT_SEQUENCE'],
                      7 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'ALT_p', 'ALT_LEN_SEQ', 'ALT_SEQUENCE']}


print(f" - {len(FaSeq_MERGE_VERSIONS)} MERGING VERSIONS IDENTIFIED :")
print("")

i = 0 

while i< len(FaSeq_MERGE_VERSIONS) :
    
    print(f"FaSeq VERSION #{i} : {', '.join(FaSeq_MERGE_VERSIONS[i])}")
    i += 1

print("")
print("")

i = 0

while i < len(FaSeq_MERGE_VERSIONS) :
    
    print(f"##### FaSeq MERGE VERSION #{i}")
    print("############################################################")
    
    # DATAFRAME CREATION BY MERGING df_REF & df_ALT
    df_FaSeq_MERGE = pd.DataFrame(df_FaSeq_REF.merge(df_FaSeq_ALT, left_on = FaSeq_REF_VERSIONS[i], right_on = FaSeq_ALT_VERSIONS[i] , how = 'outer', indicator = True))
    
    # CHANGING VALUES INSIDE 'PRESENCE' COLUMN
    df_FaSeq_MERGE._merge = ((df_FaSeq_MERGE._merge.str.replace('left_only', 'REF')).replace('right_only', 'ALT')).replace('both', 'BOTH')
    
    # COLUMN CREATION '△_PTMscore' TO VISUALIZE THE VALUES DIFFERENCES BETWEEN REF & ALT
    df_FaSeq_MERGE['p_MATCH'] = np.where(df_FaSeq_MERGE['REF_p'] == df_FaSeq_MERGE['ALT_p'], 'IDENTICAL', 'DIFFRENT')
    df_FaSeq_MERGE['LEN_SEQ_MATCH'] = np.where(df_FaSeq_MERGE["REF_LEN_SEQ"] == df_FaSeq_MERGE['ALT_LEN_SEQ'], 'IDENTICAL', 'DIFFRENT') 
    df_FaSeq_MERGE['SEQUENCE_MATCH'] = np.where(df_FaSeq_MERGE["REF_SEQUENCE"] == df_FaSeq_MERGE['ALT_SEQUENCE'], 'IDENTICAL', 'DIFFRENT') 
    
    # COLUMNS RE-ORDER
    df_FaSeq_MERGE = df_FaSeq_MERGE[['ID', 'REF_p', 'ALT_p', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'REF_LEN_SEQ', 'ALT_LEN_SEQ', 'REF_SEQUENCE', 'ALT_SEQUENCE', '_merge', 'p_MATCH', 'LEN_SEQ_MATCH', 'SEQUENCE_MATCH']]
    
    # RENAME '_merge' COLUMN BY COLUMN
    MERGE_COLUMN = f"MERGE#{i}_STATUS"
    df_FaSeq_MERGE = df_FaSeq_MERGE.rename(columns={"_merge": MERGE_COLUMN})
    
    #ADD NEW VARIABLES 
    FaSeq_BOTH = df_FaSeq_MERGE[df_FaSeq_MERGE[MERGE_COLUMN] == 'BOTH']
    FaSeq_REF = df_FaSeq_MERGE[df_FaSeq_MERGE[MERGE_COLUMN] == 'REF']
    FaSeq_ALT = df_FaSeq_MERGE[df_FaSeq_MERGE[MERGE_COLUMN] == 'ALT']
    
    OUTPUT_FaSeq_MERGE_NAME = f"{OUTPUT_FaSeq_FILE}_MERGE#{i}"
    OUTPUT_FaSeq_MERGE_FILE = f"{OUTPUT_FaSeq_MERGE_NAME}-{len(FaSeq_BOTH)}-{len(FaSeq_REF)}-{len(FaSeq_ALT)}{OUTPUT_FaSeq_FORMAT}"
    
    df_FaSeq_MERGE.to_csv(OUTPUT_FaSeq_MERGE_FILE, sep = '\t', mode='a', index=False)
    
    print("")
    print(f" - NUMBER OF SAME ROW (BOTH) : {len(FaSeq_BOTH)}")
    print(f" - NUMBER OF ROWS ONLY CONTAINED IN REF (REF) : {len(FaSeq_REF)}")
    print(f" - NUMBER OF ROWS ONLY CONTAINED IN ALT (ALT) : {len(FaSeq_ALT)}")
    print("")
    print(f" - FILE '{OUTPUT_FaSeq_MERGE_FILE}' CREATED")
    print("")
    print("")

    i += 1

print("############################################################")
print("##### PROGRAM IS DONE")
print("############################################################")
print("")