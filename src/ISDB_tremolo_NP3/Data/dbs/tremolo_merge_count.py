from __future__ import print_function

import csv
import sys

csv.field_size_limit(sys.maxsize)

if len(sys.argv) < 3:
    print(" Incorrect number of arguments")
    print(" 1. result file; 2. max results to merge and 3. and 4. count output files")
    sys.exit()

result_file = sys.argv[1]
max_results = int(sys.argv[2])

with open(result_file, 'rt') as f:
    reader = csv.reader(f, delimiter=',')
    header = False
    result_smile = {}
    result_id = {}
    result_cn = {}
    result_score = {}
    result_ppmError = {}
    result_sharedPeaks = {}
    result_mf = {}
    result_mw = {}
    result_cas = {}
    scan_count = {}
    for row in reader:
        if not header:
            scan_id_pos = row.index('msclusterID')
            unpd_id_pos = row.index('CompoundName')
            smiles_id_pos = row.index('SMILES')
            cas_id_pos = row.index('CAS')
            cn_id_pos = row.index('chemicalNames')
            mf_id_pos = row.index('molecularFormula')
            mw_id_pos = row.index('molecularWeight')
            score_id_pos = row.index('MQScore')
            ppmError_id_pos = row.index('mzErrorPPM')
            sharPeaks_id_pos = row.index('LibSearchSharedPeaks')
            header = True
        else:
            # obtain scan results
            if row[scan_id_pos] in result_smile:
                # Handle multiple molecules by node until max
                if scan_count[row[scan_id_pos]] >= max_results:
                    continue
                result_smile[row[scan_id_pos]] += "," + row[smiles_id_pos]
                result_cn[row[scan_id_pos]] += ";" + row[cn_id_pos].replace(";", ":")
                result_mf[row[scan_id_pos]] += ";" + row[mf_id_pos]
                result_mw[row[scan_id_pos]] += ";" + row[mw_id_pos]
                result_cas[row[scan_id_pos]] += ";" + row[cas_id_pos]
                result_id[row[scan_id_pos]] += ";" + row[unpd_id_pos]
                result_score[row[scan_id_pos]] += ";" + row[score_id_pos]
                result_ppmError[row[scan_id_pos]] += ";" + row[ppmError_id_pos]
                result_sharedPeaks[row[scan_id_pos]] += ";" + row[sharPeaks_id_pos]
                scan_count[row[scan_id_pos]] += 1
            else:
                result_smile[row[scan_id_pos]] = row[smiles_id_pos]
                result_cn[row[scan_id_pos]] = row[cn_id_pos].replace(";", ":")
                result_mf[row[scan_id_pos]] = row[mf_id_pos]
                result_mw[row[scan_id_pos]] = row[mw_id_pos]
                result_cas[row[scan_id_pos]] = row[cas_id_pos]
                result_id[row[scan_id_pos]] = row[unpd_id_pos]
                result_score[row[scan_id_pos]] = row[score_id_pos]
                result_ppmError[row[scan_id_pos]] = row[ppmError_id_pos]
                result_sharedPeaks[row[scan_id_pos]] = row[sharPeaks_id_pos]
                scan_count[row[scan_id_pos]] = 1

count_file = sys.argv[3]

print("  Merging to the count files result")

with open(count_file, 'rt') as f:
    reader = csv.reader(f, delimiter=',')
    header = False
    output = []
    for row in reader:
        if not header:
            cluster_index = row.index('msclusterID')
            header = row + ["UNPD_IDs_tremolo"] + ["SMILES_tremolo"] + ["chemicalNames_tremolo"] + [
                "molecularFormula_tremolo"] + ["molecularWeight_tremolo"] + ["CAS_tremolo"] + ["MQScore_tremolo"] + ["mzErrorPPM_tremolo"] + [
                         "numSharedPeaks_tremolo"]
        else:
            temp = row

            if row[cluster_index] in result_id:
                temp += [result_id[row[cluster_index]]]
            else:
                temp += ["NA"]

            if row[cluster_index] in result_smile:
                temp += [result_smile[row[cluster_index]]]
            else:
                temp += ["NA"]

            if row[cluster_index] in result_cn:
                temp += [result_cn[row[cluster_index]]]
            else:
                temp += ["NA"]

            if row[cluster_index] in result_mf:
                temp += [result_mf[row[cluster_index]]]
            else:
                temp += ["NA"]

            if row[cluster_index] in result_mw:
                temp += [result_mw[row[cluster_index]]]
            else:
                temp += ["NA"]

            if row[cluster_index] in result_cas:
                temp += [result_cas[row[cluster_index]]]
            else:
                temp += ["NA"]

            if row[cluster_index] in result_score:
                temp += [result_score[row[cluster_index]]]
            else:
                temp += ["NA"]

            if row[cluster_index] in result_ppmError:
                temp += [result_ppmError[row[cluster_index]]]
            else:
                temp += ["NA"]

            if row[cluster_index] in result_sharedPeaks:
                temp += [result_sharedPeaks[row[cluster_index]]]
            else:
                temp += ["NA"]
            output += [temp]

print("    - Outputing ", count_file)
with open(count_file, 'w') as tsvfile:
    writer = csv.writer(tsvfile, delimiter=',',
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header)
    for row in output:
        writer.writerow(row)

# also writes output for other count files (e.g. spectra count)
if len(sys.argv) == 5:
    count_file = sys.argv[4]
    print("    - Outputing ", count_file)
    with open(count_file, 'rt') as f:
        reader = csv.reader(f, delimiter=',')
        header = False
        output = []
        for row in reader:
            if not header:
                cluster_index = row.index('msclusterID')
                header = row + ["UNPD_IDs_tremolo"] + ["SMILES_tremolo"] + ["chemicalNames_tremolo"] + [
                    "molecularFormula_tremolo"] + ["molecularWeight_tremolo"] + ["CAS_tremolo"] + ["MQScore_tremolo"] + ["mzErrorPPM_tremolo"] + [
                             "numSharedPeaks_tremolo"]
            else:
                temp = row

                if row[cluster_index] in result_id:
                    temp += [result_id[row[cluster_index]]]
                else:
                    temp += ["NA"]

                if row[cluster_index] in result_smile:
                    temp += [result_smile[row[cluster_index]]]
                else:
                    temp += ["NA"]

                if row[cluster_index] in result_cn:
                    temp += [result_cn[row[cluster_index]]]
                else:
                    temp += ["NA"]

                if row[cluster_index] in result_mf:
                    temp += [result_mf[row[cluster_index]]]
                else:
                    temp += ["NA"]

                if row[cluster_index] in result_mw:
                    temp += [result_mw[row[cluster_index]]]
                else:
                    temp += ["NA"]

                if row[cluster_index] in result_cas:
                    temp += [result_cas[row[cluster_index]]]
                else:
                    temp += ["NA"]

                if row[cluster_index] in result_score:
                    temp += [result_score[row[cluster_index]]]
                else:
                    temp += ["NA"]

                if row[cluster_index] in result_ppmError:
                    temp += [result_ppmError[row[cluster_index]]]
                else:
                    temp += ["NA"]

                if row[cluster_index] in result_sharedPeaks:
                    temp += [result_sharedPeaks[row[cluster_index]]]
                else:
                    temp += ["NA"]
                output += [temp]
    with open(count_file, 'w') as tsvfile:
        writer = csv.writer(tsvfile, delimiter=',',
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header)
        for row in output:
            writer.writerow(row)
