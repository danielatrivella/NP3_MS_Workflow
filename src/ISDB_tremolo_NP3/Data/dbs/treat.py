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
    smiles_store = {}
    cn_store = {}
    mf_store = {}
    mw_store = {}
    cas_store = {}
    pm_store = {}
    for row in reader:
        if header:
            unpd_id_pos = row.index('UNPD_ID')
            smiles_pos = row.index('SMILES')
            cn_pos = row.index('cn')
            mf_pos = row.index('mf')    
            mw_pos = row.index('mw')        
            cas_pos = row.index('cas')
            pm_pos = row.index('PARENTMASS')
            header = False
        else:
            smiles_store[row[unpd_id_pos]] = row[smiles_pos]
            cn_store[row[unpd_id_pos]] = row[cn_pos].replace("\"","\'")
            mf_store[row[unpd_id_pos]] = row[mf_pos]
            mw_store[row[unpd_id_pos]] = row[mw_pos]
            cas_store[row[unpd_id_pos]] = row[cas_pos]
            pm_store[row[unpd_id_pos]] = row[pm_pos]

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
            header = row + ["SMILES"] + ["chemicalNames"] + ["molecularFormula"] + ["molecularWeight"] + ["CAS"] + ["PARENTMASS"]
        else:
            # store smiles and compoundName in the results
            temp = row   
            temp += [smiles_store[row[unpd_id_pos]]]
            temp += [cn_store[row[unpd_id_pos]]]
            temp += [mf_store[row[unpd_id_pos]]]
            temp += [mw_store[row[unpd_id_pos]]]
            temp += [cas_store[row[unpd_id_pos]]]
            temp += [pm_store[row[unpd_id_pos]]]
            output += [temp]

print("\nConverting tremolo result to CSV")       
with open(output_path+"/tremolo_results.csv", 'w') as tsvfile:
    writer = csv.writer(tsvfile, delimiter=',',
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header)
    for row in output:
        writer.writerow(row)
