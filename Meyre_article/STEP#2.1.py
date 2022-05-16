###################################################
### IMPORT PACKAGES
###################################################

import pandas as pd
import os


All_Files = ["part1.txt", "part2.txt", "part3.txt", "part4.txt", "part5.txt", "part6.txt", "part7.txt", "part8.txt", "part9.txt", "part10.txt"]

print("LIST OF FILES SELECTED TO BE MERGED :")

for File in All_Files :
    print(File)

print("")

O_File = input('INPUT THE OUTPUT_FILE : ')
print("")

df_final = pd.DataFrame()

for File in All_Files :
    
    print(File)
    print("")
    
    df_File = pd.read_csv(File, delimiter = '\t')
    
    df_File_Size = os.path.getsize(File)
    
    len_ini = len(df_File)
    
    print(f"--- MAIN INFORMATIONS ABOUT {File} ---")
    print("")
    print(f"- SIZE : {df_File_Size} Bytes ({round((df_File_Size/1000000),1)} MB)")
    print(f"- {len_ini} ROWS DETECTED")
    print("")
    print("")

    df_final = df_final.append(df_File, ignore_index = True)
    
print(df_final)

df_final.to_csv(O_File, sep="\t", mode='a', index=False) 

print("##### END #####")
