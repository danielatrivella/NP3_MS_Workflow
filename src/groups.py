#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 15:53:24 2020

@author: np3
"""

# Library importing
import pandas as pd
import argparse
import sys

#pd.set_option("precision", 4)

#---------------------------------------
# Functions
#---------------------------------------

def rows_to_skip(csv_path, sep=",", key_string='msclusterID'):
    '''
    Count how many rows to skip till find key_string parameter

    Parameters
    ----------
    csv_path : str
        Paht for count np3 output.
    sep : str, optional
        separator char of the csv. The default is ",".
    key_string : str, optional
        key word for column name identifier. The default is 'msclusterID'.

    Returns
    -------
    i : int
        Rows to skip on pandas import
        returns -1 if i == n_rows of the file.

    '''
    
    key_string = str(key_string)
    
    # No ' or "
    key_str_0 = key_string.strip("\"\'")
    # with '
    key_str_1 = '\''+key_string.strip("\"\'")+'\''
    #with "
    key_str_2 = '\"'+key_string.strip("\"\'")+'\"'
    
    f = open(csv_path).readlines()
    for i in range(len(f)):
        line = f[i].split(sep)
        if (key_str_0 in line) or \
           (key_str_1 in line) or \
           (key_str_2 in line):
               break
    
    if i == len(f)-1:
        return -1
    
    return i

#-----------------------------------------------------------------------------#

# Function to sum the group
def sum_groups(col_group, meta, count, verbose=True):

    #col_group = "GR_BAC"
    #verbose=verbose
    #meta = dfMeta
    #count = dfCount
    
    # Subgoups of GR_column without nan's
    sub_groups =  list(meta[meta[col_group].notna()][col_group].unique())
    
    if verbose:
        print("Subgroups for "+col_group)
        print(sub_groups)

    # Group Name Prefix
    pre_g_name = col_group[3:]+"_"


    for sub in sub_groups:
        
        if verbose:
            print("Computing column \""+str(pre_g_name)+str(sub)+"\"")
        
        samples = list(meta[meta[col_group]==sub][smp_cd])
        
        group_col_count = []
        for sam in samples:
            for col in count.columns:
                if (sam+"_area" == col or sam+"_spectra" == col) and \
                   (col.endswith("_area") or col.endswith("_spectra")):
                    group_col_count.append(col)
        
        # Verifying empty groups
        if len(group_col_count) > 0:
            auxCount = count[group_col_count].copy()
            count[str(pre_g_name)+str(sub)] = auxCount.sum(axis=1)
        else:
            print("Could not find columns prefixed as "+samples[0]+"_... ")
            print("Verify if some rows is needed to skip")
    
    if verbose:
        print("DONE for "+col_group, end="\n\n")

    return True

#-----------------------------------------------------------------------------#

#Function to prepend bioactivity lines
def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

#-----------------------------------------------------------------------------#
        
def arg_to_bool(string):
    string = str(string)
    if (string.lower() == 'false') or (string.lower() == '0'):
        return False
    else:
        return True


#---------------------------------------
# Parser
#---------------------------------------

# Parser for arguments
parser = argparse.ArgumentParser()

parser.add_argument('-m', '--metadata',
                    help='Path for NP3 Metadata csv file',
                    type = str, required=True)

parser.add_argument('-c', '--count_file_path',
                    help='Path for .csv clean count table without corr',
                    type = str, required=True)

parser.add_argument('-q', '--table_overwrite',
                    help="If True, overwrites the count table. Default: False",
                    type = str, default = 'False')

parser.add_argument('-v', '--verbose',
                    help="If True, show the script output information. Default: True",
                    type = str, default = 'True')



#---------------------------------------
# Parser compiling
#---------------------------------------

args = parser.parse_args()

csvMeta = args.metadata
csvCount = args.count_file_path

verbose = arg_to_bool(args.verbose)
overwrite = arg_to_bool(args.table_overwrite)

skip_rows = rows_to_skip(csvCount)

# Error handling on count table
if skip_rows == -1:
    print("\nHeader row not found on count table.\n\tAborting groups.py")
    sys.exit(0)

# Print Parameters
print("\nParameters:")
print("\tMETADATA NAME: "+csvMeta)
print("\tCOUNT TABLE NAME: "+csvCount)
print("\tSKIP ROWS ON COUNT TABLES: "+str(skip_rows))
print("\tOVERWRITE: "+str(overwrite))
print("\tVERBOSE: "+str(verbose))
print("\n")

#---------------------------------------
# Execution
#---------------------------------------


# Defining dataframes
if verbose:
    print("Creating DataFrame for process: ")

dfMeta = pd.read_csv(csvMeta)
dfCount = pd.read_csv(csvCount, skiprows=skip_rows, low_memory=False)

if verbose:
    print("DONE!")

#Verifying group col
groups_cols = []
group_pass = False
smp_cd = "SAMPLE_CODE"

if verbose:
    print("Verifying groups columns existence")

for col in dfMeta.columns:
    if "gr_"in col[:3].lower():
        groups_cols.append(col)
    
    # Get the right way to spell the sample_code column
    if "sample_code" == col.lower():
        smp_cd = col
        

if len(groups_cols) > 0:
    group_pass = True
    if verbose:
        print("Found "+str(len(groups_cols))+" group(s)")
        print(groups_cols, end="\n\n")
else:
    print("\nNo group column (GR_) found.\n\tAborting groups.py")
    sys.exit(0)


if verbose:
    print("DONE!")
       

#Function to prepend bioactivity lines
def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)
        


# Applying function
if group_pass:
    for gr in groups_cols:
        sum_groups(col_group=gr,
                   meta=dfMeta,
                   count=dfCount,
                   verbose=verbose)

    # Exporting
    if verbose:
        print("Exporting table with new groups: ")
    
    if overwrite:
        name = csvCount
    else:
        name = csvCount[:len(csvCount)-4]
        name += "_WithGroups.csv"
    dfCount.to_csv(name, index=False)
    
    #If skiped, read the correlation header
    if skip_rows > 0:
        if verbose:
            print("Preppending {} lines of bioactivity".format(skip_rows))
        #Recover lines
        append_lines = open(csvCount).readlines()[:skip_rows]
        for i in range(len(append_lines)):
            line_prepender(name, append_lines[-1-i])

    
    if verbose:
        print("DONE!")
    
    if verbose:
        print("All DONE!!")
