from Bio import SeqIO
import pandas as pd
import time 
import os 

FaSeq_DB = input("INPUT DATABASE FILE TO PARSE : ")
#FaSeq_DB = "UniprotKB_2022_01.fa"

OS = input("INPUT TYPE OF ORGANISM WANTED (Homo sapiens) : ")
#OS = "Homo sapiens"

list_FaSeq_DB = list(SeqIO.parse(FaSeq_DB, "fasta"))

FaSeq_DB_FORMAT = f".{FaSeq_DB.split('.')[-1]}"
FaSeq_DB_NAME = FaSeq_DB.replace(f'{FaSeq_DB_FORMAT}', '')

print("")
print(" - FILE NAME : ", FaSeq_DB_NAME)
print(f" - FORMAT : {FaSeq_DB_FORMAT}")
print(f" - NUMBER OF FASTA SEQUENCES DETECTED : {len(list_FaSeq_DB)}")
print(f" - ORGANISM SELETED : {OS}")
print("")

print("")
print("###########################################################################################")
print(f"##### {OS} FASTA SEQUENCES EXTRACTION FROM {FaSeq_DB}")
print("###########################################################################################")
print("")

Stime = time.time()

OUTPUT_FASTA_FORMAT = ".fa"
OUTPUT_CSV_FORMAT = ".csv"

OS_RESTRUC = "OS=" + OS + " "

OUTPUT_NAME = f"{FaSeq_DB_NAME}_{OS.replace(' ', '_')}"

FASTA_FILE = open(OUTPUT_NAME, "w")

i = 0
COUNT = 0

df_FaSeq_OS_FILE = pd.DataFrame()

while i < len(list_FaSeq_DB) : 
    
    FaSeq_DESCRIPTION = list_FaSeq_DB[i].description
    FaSeq_SEQUENCE = str(list_FaSeq_DB[i].seq)
    
    OS_ACCESSION_ID = FaSeq_DESCRIPTION.split("|")[1]
    OS_ACCESSION_NAME = (FaSeq_DESCRIPTION.split("|")[2]).split(" ")[0]
    OS_DESCRIPTION = FaSeq_DESCRIPTION.split(f"{OS_ACCESSION_NAME} ")[1]
    
    if FaSeq_DESCRIPTION.find(OS_RESTRUC) != -1 :
                      
        # df_OS_ ROW : ROW WITH ALL THE INFORMATIONS CONTAINS IN THE FASTA SEQUENCE
        df_OS_ROW = pd.DataFrame({'ACCESSION_ID' : [OS_ACCESSION_ID],
                                  'ACCESSION_NAME' : [OS_ACCESSION_NAME],
                                  'DESCRIPTION' : [OS_DESCRIPTION],
                                  'FASTA_SEQUENCE' : [FaSeq_SEQUENCE]})
        
        FASTA_FILE.write(FaSeq_DESCRIPTION + "\n")
        FASTA_FILE.write(FaSeq_SEQUENCE + "\n" + "\n")
        
        COUNT += 1
        
        df_FaSeq_OS_FILE = pd.concat([df_FaSeq_OS_FILE, df_OS_ROW])
    
    i += 1

df_FaSeq_OS_FILE = df_FaSeq_OS_FILE.reset_index(drop=True)

OUTPUT_FILE.close()    

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
print("DONE")