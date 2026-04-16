#!/usr/bin/python

import sys
import os
import pandas as pd
import ming_gnps_library
import requests
from collections import defaultdict
import argparse
import urllib.parse
#from tqdm import tqdm
from functools import lru_cache 

from tqdm.contrib.concurrent import thread_map
#import tqdm.contrib.concurrent.thread_map as thread_map
#import multiprocessing.dummy as threading

from requests.exceptions import HTTPError

# GLOBALS
# TODO: Use LRU Caching Mechanism
spectrum_id_cache = {}


@lru_cache()
def _get_gnps_library_spectrum(spectrum_id):
    # TODO: We can better sanity check the spectrum_id itself

    try:
        gnps_library_spectrum = ming_gnps_library.get_library_spectrum(spectrum_id)

        #Making sure not an error
        gnps_library_spectrum["annotations"][-1]["Compound_Name"]
    except KeyboardInterrupt:
        raise
    except HTTPError as http_err:
        print(f"Spectrum enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
        gnps_library_spectrum = None
        pass
    except:
        gnps_library_spectrum = None
        pass

    return gnps_library_spectrum

def _prep_library_dict(library_summary_df):
    library_dict = {}
    for row in library_summary_df.to_dict(orient="records"):
        library_dict[row["spectrum_id"]] = row
    return library_dict


def _enrich_gnps_annotation(output_result_dict):
    spectrum_id = output_result_dict["spectrum_id"]

    gnps_library_spectrum = _get_gnps_library_spectrum(spectrum_id)

    if gnps_library_spectrum is None:
        return output_result_dict

    # retrieving the last annotations update, not the oldest and first one in the list
    output_result_dict["Compound_Name"] = (gnps_library_spectrum["annotations"][-1]["Compound_Name"].replace("\t", ""))
    output_result_dict["Ion_Source"] = (gnps_library_spectrum["annotations"][-1]["Ion_Source"].replace("\t", ""))
    output_result_dict["Instrument"] = (gnps_library_spectrum["annotations"][-1]["Instrument"].replace("\t", ""))
    output_result_dict["Compound_Source"] = (gnps_library_spectrum["annotations"][-1]["Compound_Source"].replace("\t", ""))
    output_result_dict["PI"] = (gnps_library_spectrum["annotations"][-1]["PI"].replace("\t", ""))
    output_result_dict["Data_Collector"] = (gnps_library_spectrum["annotations"][-1]["Data_Collector"].replace("\t", ""))
    output_result_dict["Adduct"] = (gnps_library_spectrum["annotations"][-1]["Adduct"].replace("\t", ""))
    output_result_dict["Precursor_MZ"] = (gnps_library_spectrum["annotations"][-1]["Precursor_MZ"].replace("\t", ""))
    output_result_dict["ExactMass"] = (gnps_library_spectrum["annotations"][-1]["ExactMass"].replace("\t", ""))
    output_result_dict["Charge"] = (gnps_library_spectrum["annotations"][-1]["Charge"].replace("\t", ""))
    output_result_dict["CAS_Number"] = (gnps_library_spectrum["annotations"][-1]["CAS_Number"].replace("\t", ""))
    output_result_dict["Pubmed_ID"] = (gnps_library_spectrum["annotations"][-1]["Pubmed_ID"].replace("\t", ""))
    output_result_dict["Smiles"] = (gnps_library_spectrum["annotations"][-1]["Smiles"].replace("\t", ""))
    output_result_dict["INCHI"] = (gnps_library_spectrum["annotations"][-1]["INCHI"].replace("\t", ""))
    output_result_dict["INCHI_AUX"] = (gnps_library_spectrum["annotations"][-1]["INCHI_AUX"].replace("\t", ""))
    output_result_dict["Library_Class"] = (gnps_library_spectrum["annotations"][-1]["Library_Class"].replace("\t", ""))
    output_result_dict["IonMode"] = (gnps_library_spectrum["annotations"][-1]["Ion_Mode"].replace("\t", ""))

    output_result_dict["Organism"] = (gnps_library_spectrum["spectruminfo"]["library_membership"])
    # set the LibMZ equal to the MGF precursormz
    output_result_dict["LibMZ"] = (output_result_dict["precursormz"])

    if gnps_library_spectrum["annotations"][-1]["Library_Class"] == "1":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-GOLD")
        output_result_dict["LibraryQualityString"] = ("Gold")
    elif gnps_library_spectrum["annotations"][-1]["Library_Class"] == "2":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-SILVER")
        output_result_dict["LibraryQualityString"] = ("Silver")
    elif gnps_library_spectrum["annotations"][-1]["Library_Class"] == "3":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Bronze")
    elif gnps_library_spectrum["annotations"][-1]["Library_Class"] == "4":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Insilico")
    elif gnps_library_spectrum["annotations"][-1]["Library_Class"] == "5":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Insilico")
    elif gnps_library_spectrum["annotations"][-1]["Library_Class"] == "10":
        output_result_dict["UpdateWorkflowName"] = ("UPDATE-SINGLE-ANNOTATED-BRONZE")
        output_result_dict["LibraryQualityString"] = ("Challenge")
    else:
        print("Invalid Library Class", gnps_library_spectrum["annotations"][-1]["Library_Class"])

    tag_list = [ (tag["tag_desc"] + "[" + tag["tag_type"] + "]") for tag in gnps_library_spectrum["spectrum_tags"]]
    tag_string = "||".join(tag_list).replace("\t", "")

    output_result_dict["tags"] = (tag_string)
    return output_result_dict


# Here we will enrich the smiles
# changed columns: "Smiles","INCHI","molecular_formula","InChIKey","InChIKey-Planar","superclass","class","subclass","npclassifier_superclass","npclassifier_class","npclassifier_pathway","library_usi"
def _enrich_annotations(output_result_dict):
    # Calculating inchi
    if len(output_result_dict["Smiles"]) > 5 and len(output_result_dict["INCHI"]) < 5:
        try:
            inchi_url = "https://structure.gnps2.org/inchi?smiles={}".format(urllib.parse.quote_plus(output_result_dict["Smiles"]), 
                                urllib.parse.quote_plus(output_result_dict["INCHI"]))
            r = requests.get(inchi_url)
            r.raise_for_status()
            output_result_dict["INCHI"] = r.text
        except HTTPError as http_err:
            print(f"Smiles enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["INCHI"] = "N/A"
        except:
            output_result_dict["INCHI"] = "N/A"
    
    # Calculating smiles
    if len(output_result_dict["Smiles"]) < 5 and len(output_result_dict["INCHI"]) > 5:
        try:
            smiles_url = "https://structure.gnps2.org/smiles?inchi={}".format(urllib.parse.quote_plus(output_result_dict["INCHI"]), 
                                urllib.parse.quote_plus(output_result_dict["Smiles"]))
            r = requests.get(smiles_url)
            r.raise_for_status()
            output_result_dict["Smiles"] = r.text
        except HTTPError as http_err:
            print(f"Inchi enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["Smiles"] = "N/A"
        except:
            output_result_dict["Smiles"] = "N/A"

    # Calculating molecular formula
    if len(output_result_dict["Smiles"]) > 5:
        try:
            formula_url = "https://structure.gnps2.org/formula?smiles={}".format(output_result_dict["Smiles"])
            r = requests.get(formula_url)
            r.raise_for_status()
            output_result_dict["molecular_formula"] = r.text
        except HTTPError as http_err:
            print(f"Molecular formula enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["molecular_formula"] = "N/A"
        except:
            output_result_dict["molecular_formula"] = "N/A"
    else:
        output_result_dict["molecular_formula"] = "N/A"
        
    # Calculating inchi key
    if len(output_result_dict["Smiles"]) < 5 and len(output_result_dict["INCHI"]) < 5:
        output_result_dict["InChIKey"] = "N/A"
        output_result_dict["InChIKey-Planar"] = "N/A"
    else:
        try:
            inchikey_url = "https://structure.gnps2.org/inchikey?smiles={}&inchi={}".format(urllib.parse.quote_plus(output_result_dict["Smiles"]), 
                                urllib.parse.quote_plus(output_result_dict["INCHI"]))
            r = requests.get(inchikey_url)
            r.raise_for_status()
            output_result_dict["InChIKey"] = r.text
            output_result_dict["InChIKey-Planar"] = r.text.split("-")[0]
        except HTTPError as http_err:
            print(f"InChIKey enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["InChIKey"] = "N/A"
            output_result_dict["InChIKey-Planar"] = "N/A"
        except:
            output_result_dict["InChIKey"] = "N/A"
            output_result_dict["InChIKey-Planar"] = "N/A"

    # Getting Classyfire "superclass","class","subclass"
    if len(output_result_dict["InChIKey"]) > 5:
        try:
            classyfire_url = "https://classyfire.gnps2.org/entities/{}.json".format(output_result_dict["InChIKey"])
            r = requests.get(classyfire_url, timeout=0.5)
            r.raise_for_status()
            classification_json = r.json()
            output_result_dict["superclass"] = classification_json["superclass"]["name"]
            output_result_dict["class"] = classification_json["class"]["name"]
            output_result_dict["subclass"] = classification_json["subclass"]["name"]
        except HTTPError as http_err:
            print(f"Classyfire enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["superclass"] = "N/A"
            output_result_dict["class"] = "N/A"
            output_result_dict["subclass"] = "N/A"
        except:
            output_result_dict["superclass"] = "N/A"
            output_result_dict["class"] = "N/A"
            output_result_dict["subclass"] = "N/A"
    else:
        output_result_dict["superclass"] = "N/A"
        output_result_dict["class"] = "N/A"
        output_result_dict["subclass"] = "N/A"

    # Getting NP Classifier "npclassifier_superclass","npclassifier_class","npclassifier_pathway"
    if len(output_result_dict["Smiles"]) > 5:
        try:
            npclassifier_url = "https://npclassifier.gnps2.org/classify?smiles={}".format(output_result_dict["Smiles"])
            r = requests.get(npclassifier_url, timeout=10)
            r.raise_for_status()
            classification_json = r.json()

            output_result_dict["npclassifier_superclass"] = "|".join(classification_json["superclass_results"])
            output_result_dict["npclassifier_class"] = "|".join(classification_json["class_results"])
            output_result_dict["npclassifier_pathway"] = "|".join(classification_json["pathway_results"])
        except HTTPError as http_err:
            print(f"NPclassifier enrichment - HTTP error occurred: {http_err}")  # e.g. 404 Client Error
            output_result_dict["npclassifier_superclass"] = "N/A"
            output_result_dict["npclassifier_class"] = "N/A"
            output_result_dict["npclassifier_pathway"] = "N/A"
        except:
            output_result_dict["npclassifier_superclass"] = "N/A"
            output_result_dict["npclassifier_class"] = "N/A"
            output_result_dict["npclassifier_pathway"] = "N/A"
    else:
        output_result_dict["npclassifier_superclass"] = "N/A"
        output_result_dict["npclassifier_class"] = "N/A"
        output_result_dict["npclassifier_pathway"] = "N/A"

    # Adding a USI
    try:
        output_result_dict["library_usi"] = "mzspec:GNPS:GNPS-LIBRARY:{}".format(output_result_dict["SpectrumID"])
    except:
        output_result_dict["library_usi"] = "No USI"

    return output_result_dict

def enrich_spectrum_annotation_parallel(result_obj):
    # Reading existing data
    spectrum_id = result_obj["spectrum_id"]
    
    # Here we will start to write the output dictionary
    output_result_dict = result_obj
    
    # Here we are going to do the enrichment
    if "CCMSLIB" in str(spectrum_id):
        output_result_dict = _enrich_gnps_annotation(output_result_dict)
    
    return output_result_dict

def enrich_smiles_annotation_parallel(result_dict):
    # Doing further enrichment
    try:
        output_result_dict = _enrich_annotations(result_dict)
    except:
        print("SMILES enrichment failed")
        output_result_dict = result_dict
    return output_result_dict

def enrich_summary_parallell_threads(library_summary_df, output_filename, num_threads=(os.cpu_count() or 1) * 5):
    #with threading.Pool(processes=num_threads) as pool:
    #    output_list = pool.map(enrich_annotation_parallel, tqdm(library_summary_df.to_dict(orient="records")))
    # enrich parallell with tqdm interface for multiple threads
    output_list = thread_map(enrich_spectrum_annotation_parallel,
                         library_summary_df.to_dict(orient="records"),
                         max_workers=num_threads)
    output_spectrum_annotation = pd.DataFrame(output_list)
    # clear some memory space
    del output_list,library_summary_df
    # get list of unique smiles to retrieve their info
    output_unique_smiles_dict = output_spectrum_annotation.loc[(~output_spectrum_annotation.Smiles.duplicated() &
                                                                ~output_spectrum_annotation.Smiles.isna()),:].to_dict(orient="records")
    output_list = thread_map(enrich_smiles_annotation_parallel,
                             output_unique_smiles_dict,
                             max_workers=min(12,num_threads))
    output_smiles_annotation = pd.DataFrame(output_list)
    del output_unique_smiles_dict,output_list
    # merge the unique smiles annotation to the spectrum annotation
    smiles_ann_cols = ["Smiles", "INCHI", "molecular_formula", "InChIKey", "InChIKey-Planar", "superclass", "class",
                       "subclass", "npclassifier_superclass", "npclassifier_class", "npclassifier_pathway", "library_usi"]
    output_spectrum_annotation = output_spectrum_annotation.loc[:, ~output_spectrum_annotation.columns.isin(smiles_ann_cols[1:])].merge(
        output_smiles_annotation.loc[:, (smiles_ann_cols)], on="Smiles", how='left')
    # TODO Check if final size matches
    # output full annotations
    output_spectrum_annotation.to_csv(output_filename, sep="\t", index=False)

'''
from tqdm import tqdm
def enrich_summary(library_summary_df, output_filename,
                    filtertostructures=False):
    
    output_list = []
    for result_obj in tqdm(library_summary_df.to_dict(orient="records")):
        #print(result_obj)
        # Reading existing data
        spectrum_id = result_obj["spectrum_id"]

        # Here we will start to write the output dictionary
        output_result_dict = result_obj

        # Here we are going to do the enrichment
        if "CCMSLIB" in str(spectrum_id):
            output_result_dict = _enrich_gnps_annotation(output_result_dict)
     
        # Doing further enrichment
        try:
            output_result_dict = _enrich_annotations(output_result_dict)
        except:
            pass

        output_list.append(output_result_dict)

    # Here we can filter based upon the structure criteria
    if filtertostructures is True:
        # Filtering only if the length of Smiles and InchI are small
        output_list = [x for x in output_list if len(x["Smiles"]) > 5 or len(x["INCHI"]) > 5]

    pd.DataFrame(output_list).to_csv(output_filename, sep="\t", index=False)
'''

def main():
    parser = argparse.ArgumentParser(description='Pulling down GNPS identifications from online API data - needs internet connection. Enrich a library summary table with online info.')
    parser.add_argument("input_librarysummary", type=str, help="Library Summary to be enriched")
    parser.add_argument("output_filename")
    parser.add_argument("--numthreads", default=(os.cpu_count() or 1) * 5, type=int, help="Number of parallel threads to use")
    
    args = parser.parse_args()

    input_library_filename = args.input_librarysummary
    output_result_filename = args.output_filename

    # Trying to load the library summary
    if not os.path.exists(input_library_filename):
        sys.exit("The library summary file does not exists.")
    try:
        library_summary_df = pd.read_csv(input_library_filename, sep="\t", low_memory=False)
    except:
        library_summary_df = None
        sys.exit("The library summary could not be loaded, it is not a valid tsv file.")
    
    #enrich_summary(library_summary_df, output_result_filename,
    #                filtertostructures=(args.filtertostructures == "1"))
    enrich_summary_parallell_threads(library_summary_df, output_result_filename, num_threads=args.numthreads)

if __name__ == "__main__":
    main()
