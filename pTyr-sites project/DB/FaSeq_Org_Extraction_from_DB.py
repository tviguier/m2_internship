from Bio import SeqIO
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

print("")
print("###########################################################################################")
print(f"##### {ORGANISM} FASTA SEQUENCES EXTRACTION FROM {FaSeq_DB}")
print("###########################################################################################")
print("")

Stime = time.time()

OUTPUT_FORMAT = ".fa"
OUTPUT_NAME = f"FaSeq_DB_{ORGANISM.replace(' ', '_')}"
OUTPUT_FILE = f"{OUTPUT_NAME}{OUTPUT_FORMAT}"

OUTPUT = open(OUTPUT_FILE, "w")

ORGANISM_mod = "[" + ORGANISM + "]"

i = 0
COUNT = 0

while i < len(FaSeq_DB_PARSE) : 
    
    DESCRIPTION = FaSeq_DB_PARSE[i].description
    FASTA_SEQUENCE = str(FaSeq_DB_PARSE[i].seq)
    
    Homo_DESCRIPTION = (DESCRIPTION.split(ORGANISM_mod)[0]).split("\x01")[-1]
    Homo_ID = Homo_DESCRIPTION.split(" ")[0]
    Homo_ID_VERSION = Homo_ID.split(".")[-1]
    Homo_ID_NAME = (Homo_DESCRIPTION.split(";")[0]).split("RecName: Full=")[-1]
    
    if DESCRIPTION.find(ORGANISM_mod) != -1:
        
        OUTPUT.write(Homo_DESCRIPTION + "\n")
        OUTPUT.write(FASTA_SEQUENCE + "\n" + "\n")
        
        COUNT += 1
        
    i += 1

OUTPUT.close()    

NEW_OUTPUT_NAME = f"{OUTPUT_NAME}_{COUNT}"
NEW_OUTPUT_FILE = f"{NEW_OUTPUT_NAME}{OUTPUT_FORMAT}"

Etime = time.time()
Rtime = Etime - Stime

os.rename(OUTPUT_FILE, NEW_OUTPUT_FILE)

print("FASTA SEQUENCES EXTRACTION DONE")
print("")
print(f" - {COUNT} FASTA SEQUENCES EXTRACTED")
print(f" - RUN TIME : {round(Rtime, 3)}")
print(f" - '{NEW_OUTPUT_FILE}' CREATED")

print("")
print("DONE")
