'''
Script from GNPS2 NextFlow workflows
https://github.com/Wang-Bioinformatics-Lab/NextflowModules/tree/a31f16297ea5d7f996edb70897179b2e190b9778/bin/library_search
'''

from pyteomics import mgf
import argparse
import pandas as pd
import logging

def main():
    parser = argparse.ArgumentParser(description='Generate a summary of a library')
    parser.add_argument('library', help='library file')
    parser.add_argument('output_summary')
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args()

    if args.debug:
        # No practical use for now, but here for posterity
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    output_list = []

    seen_titles = set()
    seen_scans  = set()

    no_scan_err_shown   = False
    dup_scan_err_shown  = False

    # use_index=False means no requirement for title or scan number uniqueness from Pyteomics
    # To protect downstream processes, spectrum_id is checked for presence and uniqueness
    # Duplicate scan numbers are ignored as they're not included in output
    with mgf.read(args.library, use_index=False) as reader: 
        for spectrum in reader:
            # We will use the spectrum ID as the key for the dictionary
            spectrum_id = spectrum['params'].get('spectrumid', spectrum['params'].get('title', "No ID"))
            if spectrum_id == "No ID":
                raise ValueError(f"Spectrum in '{args.library}' has no 'TITLE' and no 'SPECTRUMID' attribute.")
            if spectrum_id in seen_titles:
                raise ValueError(f"Spectrum in '{args.library}' has duplicate 'TITLE' or 'SPECTRUMID': '{spectrum_id}'")
            seen_titles.add(spectrum_id)

            compound_name = spectrum['params'].get('compound_name', spectrum['params'].get("name", "No Compound"))
            smiles = spectrum['params'].get('smiles', "")
            collision_energy = spectrum['params'].get('collision_energy', 0)
            instrument = spectrum['params'].get('instrument', "")
            ion_source = spectrum['params'].get('ion_source', "")
            charge = spectrum['params'].get('charge', 0)
            adduct = spectrum['params'].get('adduct', "NA")
            
            scan = spectrum['params'].get('scans')
            if not scan and not no_scan_err_shown:
                logging.warning("An MGF entry with no scan number was found in '%s'", args.library)
                no_scan_err_shown = True
            if scan and scan in seen_scans and not dup_scan_err_shown:
                logging.warning("A duplicate scan number ('%s') was found in '%s'", str(scan), args.library)
                dup_scan_err_shown = True
            seen_scans.add(scan)

            try:
                precursormz = spectrum['params'].get('pepmass', [0])[0]
            except:
                precursormz = 0

            output_dictionary = {}
            output_dictionary["spectrum_id"] = spectrum_id
            output_dictionary["compound_name"] = compound_name
            output_dictionary["smiles"] = smiles
            output_dictionary["collision_energy"] = collision_energy
            output_dictionary["instrument"] = instrument
            output_dictionary["ion_source"] = ion_source
            output_dictionary["charge"] = charge
            output_dictionary["adduct"] = adduct
            output_dictionary["precursormz"] = precursormz
            output_dictionary["scan"] = scan
            
            output_list.append(output_dictionary)

    # creating an output df
    output_df = pd.DataFrame(output_list)
    if output_df.empty:
        # Create an empty DataFrame with the correct columns
        output_df = pd.DataFrame(columns=["spectrum_id", "compound_name", "smiles", "collision_energy", "instrument", "ion_source", "charge", "adduct", "precursormz"])
    output_df.to_csv(args.output_summary, sep="\t", index=False)

if __name__ == "__main__":
    main()