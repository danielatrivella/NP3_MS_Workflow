from __future__ import print_function

import csv
import sys

csv.field_size_limit(sys.maxsize)

if len(sys.argv) < 3:
    print(" Incorrect number of arguments")
    print(" 1. result file; 2. output path; 3. db file;")
    sys.exit()

# The tremolo file
result_file = sys.argv[1]
output_path = sys.argv[2]
db_file = sys.argv[3]

# Lazy, we load everything in memory
print(" Treating the SMILES list: ", end='')
with open(db_file, 'rt') as f:
    reader = csv.reader(f)
    header = True
    # store unpd columns
    smiles_store = {}
    cn_store = {}
    mf_store = {}
    mw_store = {}
    cas_store = {}
    pm_store = {}
    inchikey_store = {}
    # store the NPClassifier and Classyfire columns : 'NPClassifier_class', 'NPClassifier_superclass', 'NPClassifier_pathway', 'NPClassifier_isglycoside', 'ClassyFire_subclass'
    npc_superclass_store = {}
    npc_class_store = {}
    npc_pathway_store = {}
    npc_isglycoside_store = {}
    cf_subclass_store = {}
    # store the NPAtlas columns : 'NPAtlas_id', 'NPAtlas_compound_names','NPAtlas_origin_type', 'NPAtlas_genus', 'NPAtlas_origin_species'
    npa_id_store = {}
    npa_compound_names_store = {}
    npa_origin_type_store = {}
    npa_genus_store = {}
    npa_origin_species_store = {}
    for row in reader:
        if header:
            # unpd columns index
            unpd_id_pos = row.index('UNPD_ID')
            smiles_pos = row.index('Canonical_Smiles')
            cn_pos = row.index('cn')
            mf_pos = row.index('mf')    
            mw_pos = row.index('mw')        
            cas_pos = row.index('cas')
            pm_pos = row.index('PARENTMASS')
            inchikey_pos = row.index('inchik')
            # NPC columns index
            npc_superclass_pos = row.index('NPClassifier_class')
            npc_class_pos = row.index('NPClassifier_superclass')
            npc_pathway_pos = row.index('NPClassifier_pathway')
            npc_isglycoside_pos = row.index('NPClassifier_isglycoside')
            cf_subclass_pos = row.index('ClassyFire_subclass')
            # NPA columns index
            npa_id_pos = row.index('NPAtlas_id')
            npa_compound_names_pos = row.index('NPAtlas_compound_names')
            npa_origin_type_pos = row.index('NPAtlas_origin_type')
            npa_genus_pos = row.index('NPAtlas_genus')
            npa_origin_species_pos = row.index('NPAtlas_origin_species')
            header = False
        else:
            smiles_store[row[unpd_id_pos]] = row[smiles_pos]
            cn_store[row[unpd_id_pos]] = row[cn_pos].replace("\"","\'")
            mf_store[row[unpd_id_pos]] = row[mf_pos]
            mw_store[row[unpd_id_pos]] = row[mw_pos]
            cas_store[row[unpd_id_pos]] = row[cas_pos]
            pm_store[row[unpd_id_pos]] = row[pm_pos]
            inchikey_store[row[unpd_id_pos]] = row[inchikey_pos]
            npc_superclass_store[row[unpd_id_pos]] = row[npc_superclass_pos]
            npc_class_store[row[unpd_id_pos]] = row[npc_class_pos]
            npc_pathway_store[row[unpd_id_pos]] = row[npc_pathway_pos]
            npc_isglycoside_store[row[unpd_id_pos]] = row[npc_isglycoside_pos]
            cf_subclass_store[row[unpd_id_pos]] = row[cf_subclass_pos]
            npa_id_store[row[unpd_id_pos]] = row[npa_id_pos]
            npa_compound_names_store[row[unpd_id_pos]] = row[npa_compound_names_pos]
            npa_origin_type_store[row[unpd_id_pos]] = row[npa_origin_type_pos]
            npa_genus_store[row[unpd_id_pos]] = row[npa_genus_pos]
            npa_origin_species_store[row[unpd_id_pos]] = row[npa_origin_species_pos]

print("treated {} compounds".format(len(smiles_store)))

print(" Treating the result file {}".format(result_file))

with open(result_file, 'rt') as f:
    reader = csv.reader(f, delimiter='\t')
    header = False
    output = []
    for row in reader:
        if not header:
            #scan_id_pos = row.index('#Scan#')
            unpd_id_pos = row.index('CompoundName')
            #score_id_pos = row.index('MQScore')
            #ppmError_id_pos = row.index('mzErrorPPM')
            #sharPeaks_id_pos = row.index('LibSearchSharedPeaks')
            row[0] = "msclusterID"
            header = row + ["SMILES"] + ["chemicalNames"] + ["molecularFormula"] + ["molecularWeight"] + ["CAS"] + \
                     ["PARENTMASS"] + ["InChIKey"] + ["NPClassifier_superclass"] + ["NPClassifier_class"] + \
                     ["ClassyFire_subclass"] + ["NPClassifier_pathway"] + ["NPClassifier_isglycoside"] + ["NPAtlas_id"] + \
                     ["NPAtlas_compound_names"] + ["NPAtlas_origin_type"] + ["NPAtlas_genus"] + ["NPAtlas_origin_species"]
        else:
            # store smiles and compoundName in the results
            temp = row   
            temp += [smiles_store[row[unpd_id_pos]]]
            temp += [cn_store[row[unpd_id_pos]]]
            temp += [mf_store[row[unpd_id_pos]]]
            temp += [mw_store[row[unpd_id_pos]]]
            temp += [cas_store[row[unpd_id_pos]]]
            temp += [pm_store[row[unpd_id_pos]]]
            temp += [inchikey_store[row[unpd_id_pos]]]
            # npc data
            temp += [npc_superclass_store[row[unpd_id_pos]]]
            temp += [npc_class_store[row[unpd_id_pos]]]
            temp += [cf_subclass_store[row[unpd_id_pos]]]
            temp += [npc_pathway_store[row[unpd_id_pos]]]
            temp += [npc_isglycoside_store[row[unpd_id_pos]]]
            # npa data
            temp += [npa_id_store[row[unpd_id_pos]]]
            temp += [npa_compound_names_store[row[unpd_id_pos]]]
            temp += [npa_origin_type_store[row[unpd_id_pos]]]
            temp += [npa_genus_store[row[unpd_id_pos]]]
            temp += [npa_origin_species_store[row[unpd_id_pos]]]
            output += [temp]

print("\nConverting tremolo result to CSV")       
with open(output_path+"/tremolo_results.csv", 'w') as tsvfile:
    writer = csv.writer(tsvfile, delimiter=',',
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header)
    for row in output:
        writer.writerow(row)
