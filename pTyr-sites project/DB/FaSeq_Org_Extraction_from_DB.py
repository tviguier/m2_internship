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
OUTPUT_NAME = f"{ORGANISM.replace(' ', '_')}_FaSeq_DB{OUTPUT_FORMAT}"
OUTPUT_FILE = open(OUTPUT_NAME, "w")

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
        
        OUTPUT_FILE.write(Homo_DESCRIPTION + "\n")
        OUTPUT_FILE.write(FASTA_SEQUENCE + "\n" + "\n")
        
        COUNT += 1
        
    i += 1

OUTPUT_FILE.close()    

NEW_OUTPUT_NAME = f"{COUNT}_{OUTPUT_NAME}"

Etime = time.time()
Rtime = Etime - Stime

os.rename(OUTPUT_NAME, NEW_OUTPUT_NAME)

print("FASTA SEQUENCES EXTRACTION DONE")
print("")
print(f" - {COUNT} FASTA SEQUENCES EXTRACTED")
print(f" - RUN TIME : {round(Rtime, 3)}")
print(f" - '{NEW_OUTPUT_NAME}' CREATED")

print("")
print("DONE")