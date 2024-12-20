"""
Project: Diversification Rates and ADE
Author: Kristína Kocáková
Description:
Create an input file (.txt) from an .xlsx file, which is then used to generate a PyRate input file (script 2.), options for species and genus level filtering
If adenn = True then the .txt file will simply not contain the collection number, which is required for the file to be used as an ADE-NN input
"""

import os
from pandas import *

def pyrate_input(path_to_database, taxonomic_rank, path_to_output, adenn = False):
    occurrences = read_csv(path_to_database)

    # SELECT ONLY THE VALID ENTRIES (considering both age and taxonomy)

    occurrences = occurrences.loc[occurrences["age_evaluation"] == "valid"]
    occurrences = occurrences.loc[occurrences["taxonomy_validation"] == "valid"]

    #CREATE AN INPUT FILE FOR SPECIES
    if taxonomic_rank == "species":
        occurrences = occurrences.loc[occurrences["rank"] == "species"]
        occurrences = occurrences.loc[occurrences["early_interval"] != "present"]

        #Only occurrences with age resolution below 15 Myr
        occurrences = occurrences.loc[occurrences["age_range"] <= 15]

        #Select only required columns

        occurrences = occurrences[["accepted_name", "status", "max_ma", "min_ma", "collection_no"]]


        occurrences = occurrences.rename(columns={"accepted_name": "taxon_name"})

        #DROP DUPLICATES (Optional, occurrences of the same taxon from the same locality with identical ages could represent the same individual or population)
        occurrences = occurrences.drop_duplicates(keep="first")

        if adenn == True:
            occurrences = occurrences[["taxon_name", "status", "max_ma", "min_ma"]]

        #WRITE OUTPUT FILE
        occurrences.to_csv(path_to_output, sep="\t", index=False)

    #CREATE AN INPUT FILE FOR GENERA
    #INCLUDING OCCURRENCES IDENTIFIED DO A GENUS LEVEL AND THE GENUS OF OCCURRENCES IDENTIFIED TO A SPECIES LEVEL
    if taxonomic_rank == "genus":
        occurrences = occurrences.loc[(occurrences["rank"] == "species") | (occurrences["rank"] == "genus")]
        occurrences = occurrences.loc[occurrences["early_interval"] != "present"]

        #APPLY ANY ADDITIONAL FILTERS (Optional), e.g.:

        #Only occurrences with age resolution below 15 Myr
        occurrences = occurrences.loc[occurrences["age_range"] <= 15]

        #Select only required columns
        occurrences = occurrences[["genus", "genus_status", "max_ma", "min_ma", "collection_no"]]

        occurrences = occurrences.rename(columns={"genus": "taxon_name"})

        #DROP DUPLICATES (Optional, occurrences of the same taxon from the same locality with identical ages could represent the same individual or population)
        occurrences = occurrences.drop_duplicates(keep="first")

        if adenn == True:
            occurrences = occurrences.loc[["genus", "genus_status", "max_ma", "min_ma"]]

        #WRITE OUTPUT FILE
        occurrences.to_csv(path_to_output, sep="\t", index=False)



path = os.getcwd()
pyrate_input(path + "/data/occurrences.csv", "species",
             path + "/data/all_species_input.txt")
