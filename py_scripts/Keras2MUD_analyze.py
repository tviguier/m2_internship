import pandas as pd
import re
import numpy as np
import time

print("")

Keras_File = input("INPUT KERAS FILE TO ANALYZE : ")
#Keras_File = "Keras2MuD_Result_general_MuIn_ENDRESULT_ens_56001801068863.filtered.snp.txt_general_phosphorylation_Y.txt"

print("")

PGP_ID = str(input("INPUT ID OF THE PGP DATABASE : "))
#PGP_ID = str(56001801068863)

print("")

df_Keras_File = pd.read_csv(Keras_File, delimiter = '\t')
df_Keras_File.columns = ['ID', 'POSITION', 'RESIDUE', 'PTMscore']

Stime = time.time()

print(f"--- MAIN INFORMATIONS ABOUT '{Keras_File}'")
print("")
print(f" - file name : {Keras_File}")
print(f" - PGP ID : {PGP_ID}")
print(f" - Number of rows : {len(df_Keras_File)}")
print("")

print("--- START ANALYZING ---")
print("")
print(" - STEP#0 : RESTRUCTURATION ...")

df_Keras_Restruc= pd.DataFrame()

df_Keras_Restruc['Gene'] = (df_Keras_File['ID'].str.split('::', expand = True)[0]).str.split('Gene.', expand=True)[1]
df_Keras_Restruc['Type'] = np.where((df_Keras_File['ID'].str.split('__', expand = True)[1]).str.split('::', expand=True)[0] == 'REF', 'REF', 'ALT')
df_Keras_Restruc['Ensembl_ID'] = (df_Keras_File['ID'].str.split('::', expand = True)[1]).str.split('__', expand=True)[0]
df_Keras_Restruc['UniProt_ID'] = df_Keras_File['ID'].str.split('|', expand = True)[1]
df_Keras_Restruc['Name_ID'] = df_Keras_File['ID'].str.split('|', expand = True)[2]

df_Keras_Restruc['type'] = (df_Keras_File['ID'].str.split('type:', expand = True)[1]).str.split(' ', expand=True)[0]
df_Keras_Restruc['length'] = (df_Keras_File['ID'].str.split('len:', expand = True)[1]).str.split(' ', expand=True)[0]
df_Keras_Restruc['sign'] = (df_Keras_File['ID'].str.split('(', expand = True)[1]).str.split(')', expand=True)[0]
df_Keras_Restruc['score'] = (df_Keras_File['ID'].str.split('score=', expand = True)[1]).str.split(',', expand=True)[0]
df_Keras_Restruc['identity'] = df_Keras_File['ID'].str.split('|', expand = True)[3]
df_Keras_Restruc['p_value'] = (df_Keras_File['ID'].str.split('|', expand = True)[4]).str.split(',', expand=True)[0]
df_Keras_Restruc['position'] = df_Keras_File['POSITION']
df_Keras_Restruc['residue'] = df_Keras_File['RESIDUE']

df_Keras_Restruc['PTM_score'] = df_Keras_File['PTMscore']

print(" - STEP#0 : RESTRUCTURATION DONE")
print("")

### SPLITING REF & ALT ROWS ###

print(" - STEP#1 : SPLITING ROWS ...")

df_Keras_REF = df_Keras_Restruc[df_Keras_Restruc['Type'] == 'REF']
df_Keras_ALT = df_Keras_Restruc[df_Keras_Restruc['Type'] == 'ALT']

print(" - STEP#1 : SPLITING ROWS DONE")
print("")

### CREATION OF RESULT TAB ###

print(" - STEP#2 : RESULTS TAB CREATION ...")

df_Keras_PTM_01 = df_Keras_Restruc[(df_Keras_Restruc['PTM_score'] >= 0.0) & (df_Keras_Restruc['PTM_score'] <= 0.1)]
df_Keras_REF_PTM_01 = df_Keras_PTM_01[df_Keras_PTM_01['Type'] == 'REF']
df_Keras_ALT_PTM_01 = df_Keras_PTM_01[df_Keras_PTM_01['Type'] == 'ALT']

df_Keras_PTM_02 = df_Keras_Restruc[(df_Keras_Restruc['PTM_score'] > 0.1) & (df_Keras_Restruc['PTM_score'] <= 0.2)]
df_Keras_REF_PTM_02 = df_Keras_PTM_02[df_Keras_PTM_02['Type'] == 'REF']
df_Keras_ALT_PTM_02 = df_Keras_PTM_02[df_Keras_PTM_02['Type'] == 'ALT']

df_Keras_PTM_04 = df_Keras_Restruc[(df_Keras_Restruc['PTM_score'] > 0.2) & (df_Keras_Restruc['PTM_score'] <= 0.4)]
df_Keras_REF_PTM_04 = df_Keras_PTM_04[df_Keras_PTM_04['Type'] == 'REF']
df_Keras_ALT_PTM_04 = df_Keras_PTM_04[df_Keras_PTM_04['Type'] == 'ALT']

df_Keras_PTM_08 = df_Keras_Restruc[(df_Keras_Restruc['PTM_score'] > 0.4) & (df_Keras_Restruc['PTM_score'] <= 0.8)]
df_Keras_REF_PTM_08 = df_Keras_PTM_08[df_Keras_PTM_08['Type'] == 'REF']
df_Keras_ALT_PTM_08 = df_Keras_PTM_08[df_Keras_PTM_08['Type'] == 'ALT']

### INTEREST RANGE ###
df_Keras_PTM_10 = df_Keras_Restruc[(df_Keras_Restruc['PTM_score'] > 0.8) & (df_Keras_Restruc['PTM_score'] <= 1.0)]
df_Keras_REF_PTM_10 = df_Keras_PTM_10[df_Keras_PTM_10['Type'] == 'REF']
df_Keras_ALT_PTM_10 = df_Keras_PTM_10[df_Keras_PTM_10['Type'] == 'ALT']
######################

print(" - STEP#2 : RESULTS TAB CREATION DONE")
print("")
print("")

print("DATAFRAME AFTER RESTRUCTURATION \n\n", df_Keras_Restruc)
print("")

df_Keras_Results_tab_temp = pd.DataFrame({'ptyr-sites_Number' : [len(df_Keras_PTM_01),
                                                           len(df_Keras_PTM_02),
                                                           len(df_Keras_PTM_04),
                                                           len(df_Keras_PTM_08),
                                                           len(df_Keras_PTM_10),
                                                           len(df_Keras_Restruc)],
                                    
                                    'PTM_Range' : ["[0 : 0.1]", "]0.1 : 0.2]", "]0.2 : 0.4]", "]0.4 : 0.8]", "]0.8 : 1.0]", "[0 : 1.0]"]})

print("RESULTS TAB \n\n", df_Keras_Results_tab_temp)
print("")
print("")

df_Keras_Results_tab = pd.DataFrame({'REF' : [len(df_Keras_REF_PTM_01),
                                        len(df_Keras_REF_PTM_02),
                                        len(df_Keras_REF_PTM_04),
                                        len(df_Keras_REF_PTM_08),
                                        len(df_Keras_REF_PTM_10),
                                        len(df_Keras_REF)],
                               
                               'ALT' : [len(df_Keras_ALT_PTM_01),
                                        len(df_Keras_ALT_PTM_02),
                                        len(df_Keras_ALT_PTM_04),
                                        len(df_Keras_ALT_PTM_08),
                                        len(df_Keras_ALT_PTM_10),
                                        len(df_Keras_ALT)],
                               
                               'PTM_Range' : ["[0 : 0.1]", "]0.1 : 0.2]", "]0.2 : 0.4]", "]0.4 : 0.8]", "]0.8 : 1.0]", "[0 : 1.0]"]})

### FINAL RESULT TAB ###
print("RESULTS TAB WITH REF & ALT pTyr sites \n\n", df_Keras_Results_tab)
print("")
print("")

### SAVING DATAFRAMES ###

Output_prefix = f"Ouput_Keras_{PGP_ID}"
Output_format = ".csv"

Output_Keras_restruc = f"{Output_prefix}_FULL_{len(df_Keras_Restruc)}{Output_format}"
Output_Keras_REF = f"{Output_prefix}_REF_{len(df_Keras_REF)}{Output_format}"
Output_Keras_ALT = f"{Output_prefix}_ALT_{len(df_Keras_ALT)}{Output_format}"
Output_Keras_Results_tab = f"{Output_prefix}_Results_tab{Output_format}"

#Output_Keras_REF_PTM_10 = f"Ouput_keras_{PGP_ID}_REF-PTM_10-{len(df_Keras_REF_PTM_10)}.csv"
#Output_Keras_ALT_PTM_10 = f"Ouput_keras_{PGP_ID}_ALT-PTM_10-{len(df_Keras_ALT_PTM_10)}.csv"

df_Keras_Restruc.to_csv(Output_Keras_restruc, sep = '\t', mode='a', index=False)
df_Keras_REF.to_csv(Output_Keras_REF, sep = '\t', mode='a', index=False)
df_Keras_ALT.to_csv(Output_Keras_ALT, sep = '\t', mode='a', index=False)
df_Keras_Results_tab.to_csv(Output_Keras_Results_tab, sep = '\t', mode='a', index=False)

#df_Keras_REF_PTM_10.to_csv(Output_Keras_REF_PTM_10, sep = '\t', mode='a', index=False)
#df_Keras_ALT_PTM_10.to_csv(Output_Keras_ALT_PTM_10, sep = '\t', mode='a', index=False)

Etime = time.time()
Rtime = Etime - Stime

print(f" - RUN TIME : {round(Rtime, 3)}")
print("")
print(f" - '{Output_Keras_restruc}' FILE CREATED")
print(f" - '{Output_Keras_REF}' FILE CREATED")
print(f" - '{Output_Keras_ALT}' FILE CREATED")
print(f" - '{Output_Keras_Results_tab}' FILE CREATED")
#print("")
#print(f" - '{Output_Keras_REF_PTM_10}' FILE CREATED")
#print(f" - '{Output_Keras_ALT_PTM_10}' FILE CREATED")
print("")
print("")
print("END OF THE PROGRAM")
print("")
