from Bio import SeqIO
import pandas as pd
import time 
import os 

DataBase = input("INPUT FASTA DATABASE FILE TO PARSE : ")
#DataBase = "SwissProt_DB.txt"

print("")

OS = input("INPUT TYPE OF ORGANISM WANTED (Homo sapiens) : ")
#OS = "Homo sapiens"

list_DataBase = list(SeqIO.parse(DataBase, "fasta"))
DataBase_FORMAT = f".{DataBase.split('.')[-1]}"
DataBase_NAME = DataBase.replace(f'{DataBase_FORMAT}', '')

print("")
print("")
print(" - FILE NAME : ", DataBase_NAME)
print(f" - FORMAT : {DataBase_FORMAT}")
print(f" - NUMBER OF FASTA SEQUENCES DETECTED : {len(list_DataBase)}")
print(f" - ORGANISM SELETED : {OS}")
print("")
print("")
print("###########################################################################################")
print(f"##### {OS} FASTA SEQUENCES EXTRACTION FROM {DataBase}")
print("###########################################################################################")
print("")

Stime = time.time()

OUTPUT_FASTA_FORMAT = ".fa"
OUTPUT_CSV_FORMAT = ".csv"

OS_RESTRUC = "[" + OS + "]"

OUTPUT_NAME = f"{DataBase_NAME}_{OS.replace(' ', '_')}"

FASTA_FILE = open(OUTPUT_NAME, "w")

i = 0
COUNT = 0

df_FaSeq_OS_FILE = pd.DataFrame()

while i < len(list_DataBase) : 
    
    FaSeq_DESCRIPTION = list_DataBase[i].description
    FaSeq_SEQUENCE = str(list_DataBase[i].seq)
    
    OS_DESCRIPTION = FaSeq_DESCRIPTION.split(OS_RESTRUC)[0].split("\x01")[-1]
    OS_ACCESSION_ID = OS_DESCRIPTION.split(".")[0]
    OS_ACCESSION_VERSION = OS_DESCRIPTION.split(" ")[0]
    OS_ACCESSION_NAME = (OS_DESCRIPTION.split(";")[0]).split("RecName: Full=")[-1]
    
    
    if FaSeq_DESCRIPTION.find(OS_RESTRUC) != -1 :
                      
        # df_ROW : ROW WITH ALL THE INFORMATIONS CONTAINS IN THE FASTA SEQUENCE
        df_OS_ROW = pd.DataFrame({'ACCESSION_ID' : [OS_ACCESSION_ID],
                                  'ACCESSION_VERSION' : [OS_ACCESSION_VERSION],
                                  'ACCESSION_NAME' : [OS_ACCESSION_NAME],
                                  'DESCRIPTION' : [OS_DESCRIPTION],
                                  'FASTA_SEQUENCE' : [FaSeq_SEQUENCE]})
        
        FASTA_FILE.write(OS_DESCRIPTION + "\n")
        FASTA_FILE.write(FaSeq_SEQUENCE + "\n" + "\n")
        
        COUNT += 1
        
        df_FaSeq_OS_FILE = pd.concat([df_FaSeq_OS_FILE, df_OS_ROW])
    
    i += 1

df_FaSeq_OS_FILE = df_FaSeq_OS_FILE.reset_index(drop=True)

FASTA_FILE.close()    

NEW_OUTPUT_FASTA_NAME = f"{OUTPUT_NAME}-{COUNT}"
OUTPUT_FASTA_FILE = f"{NEW_OUTPUT_FASTA_NAME}{OUTPUT_FASTA_FORMAT}"

OUTPUT_CSV_NAME = f"{OUTPUT_NAME}-{len(df_FaSeq_OS_FILE)}"
OUTPUT_CSV_FILE = f"{OUTPUT_CSV_NAME}{OUTPUT_CSV_FORMAT}"

os.rename(OUTPUT_NAME, OUTPUT_FASTA_FILE)
df_FaSeq_OS_FILE.to_csv(OUTPUT_CSV_FILE, sep='\t', mode='a', index=False)

Etime = time.time()
Rtime = Etime - Stime

print("FASTA SEQUENCES EXTRACTION DONE")
print("")
print(f" - {len(df_FaSeq_OS_FILE)} FASTA SEQUENCES EXTRACTED")
print(f" - RUN TIME : {round(Rtime, 3)}")
print("")
print(f" - '{OUTPUT_FASTA_FILE}' CREATED")
print(f" - '{OUTPUT_CSV_FILE}' CREATED")
print("")
print("END")
