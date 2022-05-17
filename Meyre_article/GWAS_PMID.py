import os
import pandas as pd
from art import *
import re
import art

tprint("GWAS  PMID")
print("AUTHOR : THOMAS VIGUIER")
print("")

Folder_Path = input(">>> INPUT FOLDER PATH WITH GWAS FILES TO ANALYZE (PMID/) :")

print("")
print(f" - {len(os.listdir(Folder_Path))} FILES DETECTED IN DIRECTORY '{Folder_Path}'")
print("")

for File in os.listdir(Folder_Path) :
    print(f"   * {File}")

print("")
print("")

p_value = float(input(">>> INPUT THE p-value WANTED (5e-08) : "))

print("")
print("ANALYZE START ...")

df_Final = pd.DataFrame()
df_Final_extra = pd.DataFrame()

SNPs_count = 0

for File in os.listdir(Folder_Path) :
            
    File_Path = Folder_Path + File
    
    if Folder_Path[-1:] != '/' :
        File_Path = Folder_Path + '/' + File
    
    print("")
    print(f"- FILE '{File}' SELECTED")
    print("")
    
    df_File = pd.read_csv(File_Path, delimiter = '\t')
    df_File.fillna('None', inplace=True)

    Main_Columns = ['MAPPED_GENE', 
                    'MAPPED_TRAIT',
                    'SNPS', 
                    'STRONGEST SNP-RISK ALLELE',
                    'P-VALUE', 
                    'PUBMEDID', 
                    'STUDY ACCESSION', 
                    'FIRST AUTHOR',
                    'DATE', 
                    'DISEASE/TRAIT', 
                    'INITIAL SAMPLE SIZE',
                    'REPLICATION SAMPLE SIZE']

    df_File = df_File[Main_Columns] # TURN OFF TO SEE ALL THE COLUMNS 

    S1_Columns = ['STUDY', 
                  'ADIPONECTIN SNP AND EFFECT ALLELE', 
                  'DISCOVERY / REPLICATION POPULATION', 
                  'MAPPED GENE']

    S1_Columns_extra = ['PMID',
                        'STUDY ACCESSION',
                        'YEAR',
                        'AUTHOR',                    
                        'ADIPONECTIN SNP AND EFFECT ALLELE', 
                        'SNP',
                        'P-VALUE', 
                        'DISCOVERY / REPLICATION POPULATION', 
                        'MAPPED TRAIT',
                        'MAPPED GENE']

    df_S1 = pd.DataFrame()
    df_S1_extra = pd.DataFrame()

    i = 0
    
    while i < len(df_File) :

        ###########################
        # for S1_Columns
        
        STUDY = str(df_File['PUBMEDID'].loc[i]) + " " + str(df_File['FIRST AUTHOR'].loc[i] + " (" + str(df_File['DATE'].loc[i])[0:4] + ")" )
        ANAEA = str(df_File['STRONGEST SNP-RISK ALLELE'].loc[i])
        DRP = str(df_File['INITIAL SAMPLE SIZE'].loc[i]) + " / " + str(df_File['REPLICATION SAMPLE SIZE'].loc[i])
        MAPPED_GENE = re.split(', |,| ,| , ', str(df_File['MAPPED_GENE'].loc[i]))

        ##########################
        # for S1_Columns_extra
        
        PMID = str(df_File['PUBMEDID'].iloc[i])
        STUDY_ACCESSION = str(df_File['STUDY ACCESSION'].iloc[i])
        YEAR = str(df_File['DATE'].loc[i])[0:4]
        AUTHOR = str(df_File['FIRST AUTHOR'].loc[i])
        ANAEA = str(df_File['STRONGEST SNP-RISK ALLELE'].loc[i])
        #SNP = str((df_File['SNPS'].iloc[i]).split(':')[-1:]).replace("['", '').replace("']", '') # can't add the 'rs'
        SNP = str((df_File['SNPS'].iloc[i]))
        P_VALUE =  df_File['P-VALUE'].loc[i]
        DRP = str(df_File['INITIAL SAMPLE SIZE'].loc[i]) + " / " + str(df_File['REPLICATION SAMPLE SIZE'].loc[i])
        MAPPED_TRAIT = str(df_File['MAPPED_TRAIT'])
        MAPPED_GENE = re.split(', |,| ,| , | - ', str(df_File['MAPPED_GENE'].loc[i]))

        if P_VALUE < p_value :

            if len(MAPPED_GENE) > 1 :

                j = 0 #count number of mapped genes inside the cell

                while j < len(MAPPED_GENE) :

                    df_S1 = df_S1.append({'STUDY': STUDY, 
                                          'ADIPONECTIN SNP AND EFFECT ALLELE': ANAEA, 
                                          'DISCOVERY / REPLICATION POPULATION': DRP, 
                                          'MAPPED GENE': MAPPED_GENE[j]},ignore_index = True)

                    df_S1_extra = df_S1_extra.append({'PMID': PMID,
                                                      'STUDY ACCESSION': STUDY_ACCESSION,
                                                      'YEAR': YEAR,
                                                      'AUTHOR': AUTHOR,'ADIPONECTIN SNP AND EFFECT ALLELE': ANAEA,
                                                      'SNP': SNP,'P-VALUE': P_VALUE,
                                                      'DISCOVERY / REPLICATION POPULATION': DRP,
                                                      'MAPPED TRAIT': MAPPED_TRAIT,
                                                      'MAPPED GENE': MAPPED_GENE[j]},ignore_index = True)

                    j += 1

            else :

                df_S1 = df_S1.append({'STUDY': STUDY,
                                      'ADIPONECTIN SNP AND EFFECT ALLELE': ANAEA, 
                                      'DISCOVERY / REPLICATION POPULATION': DRP, 
                                      'MAPPED GENE': MAPPED_GENE[0]},ignore_index = True)

                df_S1_extra = df_S1_extra.append({'PMID': PMID,
                                                  'STUDY ACCESSION': STUDY_ACCESSION,
                                                  'YEAR': YEAR,
                                                  'AUTHOR': AUTHOR,'ADIPONECTIN SNP AND EFFECT ALLELE': ANAEA,
                                                  'SNP': SNP,'P-VALUE': P_VALUE,
                                                  'DISCOVERY / REPLICATION POPULATION': DRP,
                                                  'MAPPED TRAIT': MAPPED_TRAIT,
                                                  'MAPPED GENE': MAPPED_GENE[0]},ignore_index = True)

        i += 1

    PMID = str(pd.unique(df_File['PUBMEDID'])).replace('[', '').replace(']', '')
    AUTHOR = str(pd.unique(df_File['FIRST AUTHOR'])).replace("['", '').replace("']", '')
    DATE = str(pd.unique(df_File['DATE'])).replace("['", '').replace("']", '')
    STUDIES = list(pd.unique(df_File["STUDY ACCESSION"].values.ravel()))
    ASSOCIATIONS = df_File['SNPS']
    DISCOVERY_POP = str(pd.unique(df_File['INITIAL SAMPLE SIZE'])).replace("['", '').replace("']", '')
    REPLICATION_POP = str(pd.unique(df_File['REPLICATION SAMPLE SIZE'])).replace("['", '').replace("']", '')

    print(f"--- MAIN INFORMATIONS ABOUT {File} ---")
    print("")
    print(f"- PMID : {PMID}")
    print(f"- AUTHOR : {AUTHOR}")
    print(f"- DATE : {DATE}")
    print("")
    print(f"- NUMBER OF STUDIES : {len(STUDIES)}")
    print(f"- NUMBER OF ASSOCIATIONS : {len(ASSOCIATIONS)}")
    print("")
    print(f"- DISCOVERY POPULATION : {DISCOVERY_POP}")
    print(f"- REPLICATION POPULATION : {REPLICATION_POP}")
    print("")
    print(f"-- RESEARCH OF SNPS WITH P-VALUE < {p_value} ... ")

    if len(df_S1) == 0 :

        print(f"--- NO SNPS WITH P-VALUE < {p_value} DETECTED IN '{File}'")
        print("")

    else :
        
        SNPs = pd.unique(df_S1['ADIPONECTIN SNP AND EFFECT ALLELE'].values.ravel())
        
        print(f"--- {len(SNPs)} SNPS WITH P-VALUE < {p_value} DETECTED")
        print("")
        print("\n".join(SNPs))
        print("")
                     
        SNPs_count = SNPs_count + len(SNPs)

        df_Final = df_Final.append(df_S1)
        df_Final_extra = df_Final_extra.append(df_S1_extra)
    
    print("----------------------------------------------------------------------------------------")
    print("----------------------------------------------------------------------------------------")

O_File_csv = f"GWAS_SNPS_from_{Folder_Path.split('/')[0]}_p-{p_value}.csv"
O_File_extra_csv = f"GWAS_SNPS_extra_from_{Folder_Path.split('/')[0]}_p-{p_value}.csv"

df_Final.to_csv(O_File_csv, sep="\t", mode='a', index=False)
df_Final_extra.to_csv(O_File_extra_csv, sep="\t", mode='a', index=False)

LDProxy_SNPS = list(pd.unique(df_Final_extra["SNP"].values.ravel()))

O_File_txt = "GWAS_list_SNPS.txt"
O_File = open(O_File_txt, 'w+')

O_File.write(f"All SNPS from : {Folder_Path} \n" + "\n")
O_File.write("\t" + "\n \t".join(os.listdir(Folder_Path)) + "\n" + "\n")
O_File.write("All these SNPS can be parse with LDProxy algorithm to visualize their respective proxy SNPS. \n" + "\n")
O_File.write("\n".join(LDProxy_SNPS))

O_File.close()

print("")
print("")
print(f"TOTAL OF {SNPs_count} UNIQUE SNPS WITH P-VALUE < {p_value} DETECTED")
print("")
print(f"- FILE '{O_File_csv}' AND '{O_File_extra_csv}' WERE CREATED WITH ALL THE SNPS")
print(f"- FILE '{O_File_txt}' WAS CREATED WITH ALL THE SNPS ID TO BE PARSE IN LDPROX ALGORITHM")
print("")
print("##### END #####")
