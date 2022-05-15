"""
pTyr RESTRUCTURATION
AUTHOR : THOMAS VIGUIER
AIM : 
- CONVERT TAB.txt FILE WITH pTyr ROWS TO DATAFRAME
- SPLITING REF & ALT ROWS INTO 2 DATAFRAMES
- MERGE THE DATAFRAMES TO SEE THE COMMON AND DIFFERENT pTyr SITES
"""

import pandas as pd
import numpy as np
import re
import time


POS_INPUTS = ['YES', 'Y', 'yes', 'y', 'oui', 'OUI']

print("")
print("########################################################################################################################")
print("########################################################################################################################")
print("")
print("############### pTyr SITES ANALYZE")
print("")
print("########################################################################################################################")
print("########################################################################################################################")
print("")

pTyr_FILE = input(">>> INPUT pTyr FILE : ") 

print("")

COMMON_ID = input(f">>> INPUT ID OF '{pTyr_FILE}' : ")

print("")
print("###########################################################################################")
print("###########################################################################################")
print("")
print("########## I : CREATION OF pTyr-REF & pTyr-ALT FILES")
print("")
print("###########################################################################################")
print("###########################################################################################")
print("")
    
print(f"##### FILE '{pTyr_FILE}' SELECTED")

print("")
print("###########################################################################################")
print(f"##### STEP I.1 : MAIN INFORMATIONS ABOUT {pTyr_FILE}")
print("###########################################################################################")
print("")

df_pTyr_FILE = pd.read_csv(pTyr_FILE, delimiter = '\t')
pTyr_FORMAT = f".{pTyr_FILE.split('.')[-1]}"
pTyr_FILE_NAME = pTyr_FILE.replace(f'.{pTyr_FORMAT}', '')

print(" - FILE NAME : ", pTyr_FILE_NAME)
print(" - ID : ", COMMON_ID)
print(f" - FORMAT : {pTyr_FORMAT}")
print(f" - NUMBER OF pTyr-Sites DETECTED : {len(df_pTyr_FILE)}")
print("")

df_pTyr_FILE = pd.read_csv(pTyr_FILE, delimiter = '\t')
df_pTyr_FILE.columns =['ID', 'POSITION', 'RESIDUE', 'PTMscore']

OUTPUT_pTyr_FORMAT = ".csv" 

OUTPUT_pTyr_NAME = f"pTyr_{COMMON_ID}-{len(df_pTyr_FILE)}"
OUTPUT_pTyr_FILE = f"{OUTPUT_pTyr_NAME}{OUTPUT_pTyr_FORMAT}"

print("")
print("###########################################################################################")
print(f"##### STEP I.2 : RESTRUCTURATION")
print("###########################################################################################")
print("")

Stime = time.time()

print("RESTRUCTURATION ... (TXT TO DATAFRAME)")
print("")


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

Etime = time.time()
Rtime = Etime - Stime

print("RESTRUCTURATION DONE")
print("")
print(f" - RUN TIME : {round(Rtime, 3)}")
print("")

pTyr_INPUT_0 = input(">>> PRINT THE RESULT DATAFRAME ? (Y/N) : ")

print("")

if pTyr_INPUT_0 in POS_INPUTS :
    
    print("")
    print(df_pTyr_RESTRUC)
    print("")

pTyr_INPUT_1 = input(">>> GROUP THE COLUMNS ['TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN'] IN A NEW DATAFRAME ? (Y/N) : ")

if pTyr_INPUT_1 in POS_INPUTS :
    
    df_pTyr_HEADER = pd.DataFrame()
    
    df_pTyr_HEADER['ID'] = df_pTyr_RESTRUC['ID']
    df_pTyr_HEADER['p'] = df_pTyr_RESTRUC['p']
    df_pTyr_HEADER['HEADER'] = "type:" + df_pTyr_RESTRUC['TYPE'] + " len:" + df_pTyr_RESTRUC['LENGTH'].astype(str) + " gc:" + df_pTyr_RESTRUC['GC'] + " range:" + df_pTyr_RESTRUC['RANGE'].astype(str) + "(" + df_pTyr_RESTRUC['SIGN'] + ")"
    df_pTyr_HEADER['POSITION'] = df_pTyr_RESTRUC['POSITION']
    df_pTyr_HEADER['RESIDUE'] = df_pTyr_RESTRUC['RESIDUE']
    df_pTyr_HEADER['PTMscore'] = df_pTyr_RESTRUC['PTMscore']
    
    df_pTyr_HEADER = df_pTyr_HEADER[['ID', 'p', 'HEADER', 'POSITION', 'RESIDUE', 'PTMscore']]
    
    print("")
    print("GROUP DONE")
    print("")
    
    pTyr_INPUT_2 = input(">>> PRINT THE HEADER DATAFRAME ? (Y/N) : ")
    
    if pTyr_INPUT_2 in POS_INPUTS :
        
        print("")
        print(df_pTyr_HEADER)

print("")
print("###########################################################################################")
print(f"##### STEP I.3 : SPLITING pTyr-REF & pTyr-ALT ROWS")
print("###########################################################################################")
print("")

print("SPLITING ...")
print("")

Stime = time.time()

df_pTyr_REF = df_pTyr_RESTRUC[df_pTyr_RESTRUC['pTyr_TYPE'] == 'REF']
df_pTyr_REF = df_pTyr_REF.rename(columns={'p' : 'REF_p', 'PTMscore': 'REF_PTMscore'})
df_pTyr_REF = df_pTyr_REF[['ID', 'pTyr_TYPE', 'REF_p', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'REF_PTMscore']]
                 
df_pTyr_ALT = df_pTyr_RESTRUC[df_pTyr_RESTRUC['pTyr_TYPE'] == 'ALT']
df_pTyr_ALT = df_pTyr_ALT.rename(columns={'p' : 'ALT_p', 'PTMscore': 'ALT_PTMscore'})
df_pTyr_ALT = df_pTyr_ALT[['ID', 'pTyr_TYPE', 'ALT_p', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'ALT_PTMscore']]

Etime = time.time()
Rtime = Etime - Stime


print("SPLITING DONE")
print("")
print(f" - RUN TIME : {round(Rtime, 3)}")
print("")

pTyr_INPUT_3 = input(">>> PRINT pTyr-REF & pTyr-ALT DATAFRAMES ? (Y/N) : ")

print("")

if pTyr_INPUT_3 in POS_INPUTS :
    
    print(" - df_pTyr_REF\n", df_pTyr_REF)
    print("")
    print(" - df_pTyr_ALT\n",df_pTyr_ALT)

print("")

pTyr_INPUT_4 = input(">>> SAVE pTyr_RESULTS DATAFRAMES ? (Y/N) : ")

print("")

OUTPUT_pTyr_REF_NAME = f"{OUTPUT_pTyr_NAME}_REF-{len(df_pTyr_REF)}"
OUTPUT_pTyr_REF_FILE = f"{OUTPUT_pTyr_REF_NAME}{OUTPUT_pTyr_FORMAT}"

OUTPUT_pTyr_ALT_NAME = f"{OUTPUT_pTyr_NAME}_ALT-{len(df_pTyr_ALT)}"
OUTPUT_pTyr_ALT_FILE = f"{OUTPUT_pTyr_ALT_NAME}{OUTPUT_pTyr_FORMAT}"

if pTyr_INPUT_4 in POS_INPUTS :

    df_pTyr_FILE.to_csv(f"{OUTPUT_pTyr_FILE}", sep='\t', mode='a', index=False)
    df_pTyr_REF.to_csv(OUTPUT_pTyr_REF_FILE, sep='\t', mode='a', index=False)
    df_pTyr_ALT.to_csv(OUTPUT_pTyr_ALT_FILE, sep='\t', mode='a', index=False)
    
    print(f" - FILE '{OUTPUT_pTyr_FILE}' CREATED")
    
    if pTyr_INPUT_1 in POS_INPUTS :
        
        OUTPUT_pTyr_HEADER_NAME = f"{OUTPUT_pTyr_NAME}_HEADER"
        OUTPUT_pTyr_HEADER_FILE = f"{OUTPUT_pTyr_HEADER_NAME}{OUTPUT_pTyr_FORMAT}"
        
        df_pTyr_HEADER.to_csv(OUTPUT_pTyr_HEADER_FILE, sep='\t', mode='a', index=False)
        
        print(f" - FILE '{OUTPUT_pTyr_HEADER_FILE}' CREATED")
          
    print(f" - FILE '{OUTPUT_pTyr_REF_FILE}' CREATED")
    print(f" - FILE '{OUTPUT_pTyr_ALT_FILE}' CREATED")
    print("")


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

print("")
print("")
print("###########################################################################################")
print(f"##### VIJAY TAB")
print("###########################################################################################")
print("")

df_pTyr_REF = pd.read_csv(OUTPUT_pTyr_REF_FILE, delimiter = '\t')
df_pTyr_ALT = pd.read_csv(OUTPUT_pTyr_ALT_FILE, delimiter = '\t')

Stime = time.time()

df_pTyr_REF_01 = pd.DataFrame(df_pTyr_REF[(df_pTyr_REF['REF_PTMscore'] >= 0.0) & (df_pTyr_REF['REF_PTMscore'] <= 0.1)])
df_pTyr_ALT_01 = pd.DataFrame(df_pTyr_ALT[(df_pTyr_ALT['ALT_PTMscore'] >= 0.0) & (df_pTyr_ALT['ALT_PTMscore'] <= 0.1)])

OUTPUT_pTyr_REF_01_NAME = f"{OUTPUT_pTyr_REF_NAME}_01-{len(df_pTyr_REF_01)}"
OUTPUT_pTyr_REF_01_FILE = f"{OUTPUT_pTyr_REF_01_NAME}{OUTPUT_pTyr_FORMAT}" 
OUTPUT_pTyr_ALT_01_NAME = f"{OUTPUT_pTyr_ALT_NAME}_01-{len(df_pTyr_ALT_01)}"
OUTPUT_pTyr_ALT_01_FILE = f"{OUTPUT_pTyr_ALT_01_NAME}{OUTPUT_pTyr_FORMAT}" 
df_pTyr_REF_01.to_csv(OUTPUT_pTyr_REF_01_FILE, sep='\t', mode='a', index=False)
df_pTyr_ALT_01.to_csv(OUTPUT_pTyr_ALT_01_FILE, sep='\t', mode='a', index=False)

###############

df_pTyr_REF_02 = pd.DataFrame(df_pTyr_REF[(df_pTyr_REF['REF_PTMscore'] > 0.1) & (df_pTyr_REF['REF_PTMscore'] <= 0.2)])
df_pTyr_ALT_02 = pd.DataFrame(df_pTyr_ALT[(df_pTyr_ALT['ALT_PTMscore'] > 0.1) & (df_pTyr_ALT['ALT_PTMscore'] <= 0.2)])

OUTPUT_pTyr_REF_02_NAME = f"{OUTPUT_pTyr_REF_NAME}_02-{len(df_pTyr_REF_02)}"
OUTPUT_pTyr_REF_02_FILE = f"{OUTPUT_pTyr_REF_02_NAME}{OUTPUT_pTyr_FORMAT}" 
OUTPUT_pTyr_ALT_02_NAME = f"{OUTPUT_pTyr_ALT_NAME}_02-{len(df_pTyr_ALT_02)}"
OUTPUT_pTyr_ALT_02_FILE = f"{OUTPUT_pTyr_ALT_02_NAME}{OUTPUT_pTyr_FORMAT}" 
df_pTyr_REF_02.to_csv(OUTPUT_pTyr_REF_02_FILE, sep='\t', mode='a', index=False)
df_pTyr_ALT_02.to_csv(OUTPUT_pTyr_ALT_02_FILE, sep='\t', mode='a', index=False)

###############

df_pTyr_REF_04 = pd.DataFrame(df_pTyr_REF[(df_pTyr_REF['REF_PTMscore'] > 0.2) & (df_pTyr_REF['REF_PTMscore'] <= 0.4)])
df_pTyr_ALT_04 = pd.DataFrame(df_pTyr_ALT[(df_pTyr_ALT['ALT_PTMscore'] > 0.2) & (df_pTyr_ALT['ALT_PTMscore'] <= 0.4)])

OUTPUT_pTyr_REF_04_NAME = f"{OUTPUT_pTyr_REF_NAME}_04-{len(df_pTyr_REF_04)}"
OUTPUT_pTyr_REF_04_FILE = f"{OUTPUT_pTyr_REF_04_NAME}{OUTPUT_pTyr_FORMAT}" 
OUTPUT_pTyr_ALT_04_NAME = f"{OUTPUT_pTyr_ALT_NAME}_04-{len(df_pTyr_ALT_04)}"
OUTPUT_pTyr_ALT_04_FILE = f"{OUTPUT_pTyr_ALT_04_NAME}{OUTPUT_pTyr_FORMAT}" 
df_pTyr_REF_04.to_csv(OUTPUT_pTyr_REF_04_FILE, sep='\t', mode='a', index=False)
df_pTyr_ALT_04.to_csv(OUTPUT_pTyr_ALT_04_FILE, sep='\t', mode='a', index=False)

###############

df_pTyr_REF_08 = pd.DataFrame(df_pTyr_REF[(df_pTyr_REF['REF_PTMscore'] > 0.4) & (df_pTyr_REF['REF_PTMscore'] <= 0.8)])
df_pTyr_ALT_08 = pd.DataFrame(df_pTyr_ALT[(df_pTyr_ALT['ALT_PTMscore'] > 0.4) & (df_pTyr_ALT['ALT_PTMscore'] <= 0.8)])

OUTPUT_pTyr_REF_08_NAME = f"{OUTPUT_pTyr_REF_NAME}_08-{len(df_pTyr_REF_08)}"
OUTPUT_pTyr_REF_08_FILE = f"{OUTPUT_pTyr_REF_08_NAME}{OUTPUT_pTyr_FORMAT}" 
OUTPUT_pTyr_ALT_08_NAME = f"{OUTPUT_pTyr_ALT_NAME}_08-{len(df_pTyr_ALT_08)}"
OUTPUT_pTyr_ALT_08_FILE = f"{OUTPUT_pTyr_ALT_08_NAME}{OUTPUT_pTyr_FORMAT}" 
df_pTyr_REF_08.to_csv(OUTPUT_pTyr_REF_08_FILE, sep='\t', mode='a', index=False)
df_pTyr_ALT_08.to_csv(OUTPUT_pTyr_ALT_08_FILE, sep='\t', mode='a', index=False)

###############

df_pTyr_REF_10 = pd.DataFrame(df_pTyr_REF[(df_pTyr_REF['REF_PTMscore'] > 0.8) & (df_pTyr_REF['REF_PTMscore'] <= 1.0)])
df_pTyr_ALT_10 = pd.DataFrame(df_pTyr_ALT[(df_pTyr_ALT['ALT_PTMscore'] > 0.8) & (df_pTyr_ALT['ALT_PTMscore'] <= 1.0)])

OUTPUT_pTyr_REF_10_NAME = f"{OUTPUT_pTyr_REF_NAME}_10-{len(df_pTyr_REF_10)}"
OUTPUT_pTyr_REF_10_FILE = f"{OUTPUT_pTyr_REF_10_NAME}{OUTPUT_pTyr_FORMAT}" 
OUTPUT_pTyr_ALT_10_NAME = f"{OUTPUT_pTyr_ALT_NAME}_10-{len(df_pTyr_ALT_10)}"
OUTPUT_pTyr_ALT_10_FILE = f"{OUTPUT_pTyr_ALT_10_NAME}{OUTPUT_pTyr_FORMAT}" 
df_pTyr_REF_10.to_csv(OUTPUT_pTyr_REF_10_FILE, sep='\t', mode='a', index=False)
df_pTyr_ALT_10.to_csv(OUTPUT_pTyr_ALT_10_FILE, sep='\t', mode='a', index=False)

###############

TOTAL_pTyr_REF = len(df_pTyr_REF_01 + df_pTyr_REF_02 + df_pTyr_REF_04 + df_pTyr_REF_08 + df_pTyr_REF_10)
TOTAL_pTyr_ALT = len(df_pTyr_ALT_01 + df_pTyr_ALT_02 + df_pTyr_ALT_04 + df_pTyr_ALT_08 + df_pTyr_ALT_10)

#df_pTyr_GLOBAL = pd.DataFrame()
df_pTyr_GLOBAL = pd.DataFrame({'REF': [len(df_pTyr_REF_01), len(df_pTyr_REF_02), len(df_pTyr_REF_04), len(df_pTyr_REF_08), len(df_pTyr_REF_10), TOTAL_pTyr_REF], 
                               'ALT': [len(df_pTyr_ALT_01), len(df_pTyr_ALT_02), len(df_pTyr_ALT_04), len(df_pTyr_ALT_08), len(df_pTyr_ALT_10), TOTAL_pTyr_ALT], 
                               'PTMscore': ["[0 : 0.1]", "]0.1 : 0.2]", "]0.2 : 0.4]", "]0.4 : 0.8]", "]0.8 : 1.0]", "[0 : 1.0]"]})


OUTPUT_pTyr_GLOBAL_NAME = f"{OUTPUT_pTyr_NAME}_PTMscores"
OUTPUT_pTyr_GLOBAL_FILE = f"{OUTPUT_pTyr_GLOBAL_NAME}{OUTPUT_pTyr_FORMAT}"
df_pTyr_GLOBAL.to_csv(OUTPUT_pTyr_GLOBAL_FILE, sep='\t', mode='a', index=False)

Etime = time.time()
Rtime = Etime - Stime

print(f" - REF NUMBER pTyr SITES : {len(df_pTyr_REF)}")
print(f" - ALT NUMBER pTyr SITES : {len(df_pTyr_ALT)}")
print("")
print(df_pTyr_GLOBAL)
print("")
print(f" - RUN TIME : {round(Rtime, 3)}")
print(" - ALL THE DATAFRAMES ARE SAVED")
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
print("########## II : MERGING pTyr-REF & pTyr-ALT DATAFRAMES")
print("")
print("###########################################################################################")
print("###########################################################################################")
print("")

df_pTyr_REF = pd.read_csv(OUTPUT_pTyr_REF_FILE, delimiter = '\t')
df_pTyr_ALT = pd.read_csv(OUTPUT_pTyr_ALT_FILE, delimiter = '\t')


Stime = time.time()


pTyr_MERGE_VERSIONS = {0 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE'], # COMMON MERGE
                       1 : ['ID', 'REF_p', 'ALT_p', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE'], # MERGE WITH ADDING RESPECTIVE p
                       2 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'REF_PTMscore', 'ALT_PTMscore'], # MERGE WITH ADDING RESPECTIVE PTMscore
                       3 : ['ID', 'REF_p', 'ALT_p', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'REF_PTMscore', 'ALT_PTMscore']} # MERGE WITH ADDING RESPECTIVE p & PTMscore

pTyr_REF_VERSIONS = {0 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE'],
                     1 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'REF_p'],
                     2 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'REF_PTMscore'],
                     3 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'REF_p', 'REF_PTMscore']}

pTyr_ALT_VERSIONS = {0 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE'],
                     1 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'ALT_p'],
                     2 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'ALT_PTMscore'],
                     3 : ['ID', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'ALT_p', 'ALT_PTMscore']}


print(f" - {len(pTyr_MERGE_VERSIONS)} MERGING VERSIONS IDENTIFIED :")
print("")

i = 0 

while i< len(pTyr_MERGE_VERSIONS) :
    
    print(f"pTyr VERSION #{i} : {', '.join(pTyr_MERGE_VERSIONS[i])}")
    i += 1

print("")
print("")

i = 0

while i < len(pTyr_MERGE_VERSIONS) :
    
    print(f"##### pTyr MERGE VERSION #{i}")
    print("############################################################")
    
    # DATAFRAME CREATION BY MERGING df_REF & df_ALT
    df_pTyr_MERGE = pd.DataFrame(df_pTyr_REF.merge(df_pTyr_ALT, left_on = pTyr_REF_VERSIONS[i], right_on = pTyr_ALT_VERSIONS[i] , how = 'outer', indicator = True))
    
    # CHANGING VALUES INSIDE 'PRESENCE' COLUMN
    df_pTyr_MERGE._merge = ((df_pTyr_MERGE._merge.str.replace('left_only', 'REF')).replace('right_only', 'ALT')).replace('both', 'BOTH')
    
    # COLUMN CREATION '△_PTMscore' TO VISUALIZE THE VALUES DIFFERENCES BETWEEN REF & ALT
    df_pTyr_MERGE['p_MATCH'] = np.where(df_pTyr_MERGE['REF_p'] == df_pTyr_MERGE['ALT_p'], 'IDENTICAL', 'DIFFRENT')
    df_pTyr_MERGE['△_PTMscore'] = np.where(df_pTyr_MERGE["REF_PTMscore"] == df_pTyr_MERGE['ALT_PTMscore'], 0, df_pTyr_MERGE["REF_PTMscore"] - df_pTyr_MERGE["ALT_PTMscore"]) 
    
    # COLUMNS RE-ORDER
    df_pTyr_MERGE = df_pTyr_MERGE[['ID', 'REF_p', 'ALT_p', 'TYPE', 'LENGTH', 'GC', 'RANGE', 'SIGN', 'POSITION', 'RESIDUE', 'REF_PTMscore', 'ALT_PTMscore', '_merge', 'p_MATCH', '△_PTMscore']]
    
    # RENAME '_merge' COLUMN BY COLUMN
    MERGE_COLUMN = f"MERGE#{i}_STATUS"
    df_pTyr_MERGE = df_pTyr_MERGE.rename(columns={"_merge": MERGE_COLUMN})
    
    #ADD NEW VARIABLES 
    pTyr_BOTH = df_pTyr_MERGE[df_pTyr_MERGE[MERGE_COLUMN] == 'BOTH']
    pTyr_REF = df_pTyr_MERGE[df_pTyr_MERGE[MERGE_COLUMN] == 'REF']
    pTyr_ALT = df_pTyr_MERGE[df_pTyr_MERGE[MERGE_COLUMN] == 'ALT']
    
    OUTPUT_pTyr_MERGE_NAME = f"{OUTPUT_pTyr_FILE}_MERGE#{i}"
    OUTPUT_pTyr_MERGE_FILE = f"{OUTPUT_pTyr_MERGE_NAME}-{len(pTyr_BOTH)}-{len(pTyr_REF)}-{len(pTyr_ALT)}{pTyr_FORMAT}"
    
    df_pTyr_MERGE.to_csv(OUTPUT_pTyr_MERGE_FILE, sep = '\t', mode='a', index=False)
    
    print("")
    print(f"##### MERGING#{i} RESULTS")
    print("")
    print(f" - NUMBER OF SAME ROW (BOTH) : {len(pTyr_BOTH)}")
    print(f" - NUMBER OF ROWS ONLY CONTAINED IN REF (REF) : {len(pTyr_REF)}")
    print(f" - NUMBER OF ROWS ONLY CONTAINED IN ALT (ALT) : {len(pTyr_ALT)}")
    print("")
    print(f" - FILE '{OUTPUT_pTyr_MERGE_FILE}' CREATED")
    print("")
    print("")

    i += 1

print("############################################################")
print("##### PROGRAM IS DONE")
print("############################################################")
print("")