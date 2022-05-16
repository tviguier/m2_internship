from Bio import SeqIO
import pandas as pd
import time 
import os 

FaSeq_DB = input("INPUT DATABASE FILE TO PARSE : ")
ORGANISM = input("INPUT TYPE OF ORGANISM WANTED (Homo sapiens) : ")

print("")
print("###########################################################################################")
print(f"##### MAIN INFORMATIONS ABOUT {FaSeq_DB}")
print("###########################################################################################")
print("")

FaSeq_DB_PARSE = list(SeqIO.parse(FaSeq_DB, "fasta"))

FaSeq_DB_FORMAT = f".{FaSeq_DB.split('.')[-1]}"
FaSeq_DB_NAME = FaSeq_DB.replace(f'.{FaSeq_DB_FORMAT}', '')

print(" - FILE NAME : ", FaSeq_DB_NAME)
print(f" - FORMAT : {FaSeq_DB_FORMAT}")
print(f" - NUMBER OF FASTA SEQUENCES DETECTED : {len(FaSeq_DB_PARSE)}")
print(f" - ORGANISM SELETED : {ORGANISM}")

print("")
print("###########################################################################################")
print(f"##### {ORGANISM} FASTA SEQUENCES EXTRACTION FROM {FaSeq_DB}")
print("###########################################################################################")
print("")

Stime = time.time()

OUTPUT_FASTA_FORMAT = ".fa"
OUTPUT_CSV_FORMAT = ".csv"

OUTPUT_FASTA_NAME = f"FaSeq_DB_{ORGANISM.replace(' ', '_')}"
OUTPUT_FASTA_FILE = f"{OUTPUT_NAME}{OUTPUT_FASTA_FORMAT}"

OUTPUT = open(OUTPUT_FILE, "w")

ORGANISM_mod = "[" + ORGANISM + "]"

df_DB_FILE = pd.DataFrame()

i = 0
COUNT = 0

while i < len(FaSeq_DB_PARSE) : 
    
    DESCRIPTION = str(FaSeq_DB_PARSE[i].description)
    FASTA_SEQUENCE = str(FaSeq_DB_PARSE[i].seq)
    
    Homo_DESCRIPTION = (DESCRIPTION.split(ORGANISM_mod)[0]).split("\x01")[-1]
    Homo_ID = Homo_DESCRIPTION.split(" ")[0]
    Homo_ID_VERSION = int(Homo_ID.split(".")[-1])
    Homo_ID_NAME = (Homo_DESCRIPTION.split(";")[0]).split("RecName: Full=")[-1]
    
    if DESCRIPTION.find(ORGANISM_mod) != -1:
        
        OUTPUT.write(Homo_DESCRIPTION + "\n")
        OUTPUT.write(FASTA_SEQUENCE + "\n" + "\n")
        
        df_ROW = pd.DataFrame({'ID' : [Homo_ID],
                               'ID_VERSION' : [Homo_ID_VERSION],
                               'ID_NAME' : [Homo_ID_NAME],
                               'DESCRIPTION' : [Homo_DESCRIPTION],
                               'LEN_SEQ' : [len(FASTA_SEQUENCE)],
                               'SEQUENCE' : [FASTA_SEQUENCE]})
        
        
        df_DB_FILE = pd.concat([df_DB_FILE, df_ROW])
           
        COUNT += 1
        
    i += 1

df_DB_FILE = df_DB_FILE.reset_index(drop=True)    

OUTPUT.close()    

NEW_OUTPUT_FASTA_NAME = f"{OUTPUT_FASTA_NAME}_{COUNT}"
NEW_OUTPUT_FASTA_FILE = f"{NEW_OUTPUT_FASTA_NAME}{OUTPUT_FASTA_FORMAT}"

OUTPUT_CSV_FILE = f"{NEW_OUTPUT_FASTA_NAME}{OUTPUT_CSV_FORMAT}"

df_DB_FILE.to_csv(CSV_OUTPUT_FILE, sep='\t', mode='a', index=False)

Etime = time.time()
Rtime = Etime - Stime

os.rename(OUTPUT_FASTA_FILE, NEW_OUTPUT_FASTA_FILE)

print("FASTA SEQUENCES EXTRACTION DONE")
print("")
print(f" - {COUNT} FASTA SEQUENCES EXTRACTED")
print(f" - RUN TIME : {round(Rtime, 3)}")
print(f" - '{NEW_OUTPUT_FASTA_FILE}' CREATED")
print(f" - '{OUTPUT_CSV_FILE}' CREATED")
print("")
print("DONE")
