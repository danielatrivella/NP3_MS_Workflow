# @Crisfbazz
# Adapted from https://github.com/Wang-Bioinformatics-Lab/LibrarySearch_Workflow/blob/master/data/get_data.sh
# commit 4284470
# Following the flow present in the GNPS2 LibrarySearch_Workflow
# adapted here to be executed offline and using only the ALL_GNPS_NO_PROPOGATED.mgf library file
cd src/GNPS_LibrarySearch/data

# mkdir if not created yet
mkdir -p libraries && cd libraries
# download the library file if a new version is present in the server
#wget -N https://external.gnps2.org/gnpslibrary/GNPS-LIBRARY.mgf
LIBRARY_FILE=ALL_GNPS_NO_PROPOGATED_SPLITS.mgf.tar.gz
if [ -f "$LIBRARY_FILE" ]; then
    t1=$(stat -c %y $LIBRARY_FILE)
else
    t1="no library file is present"
fi
if wget -N https://external.gnps2.org/gnpslibrary/ALL_GNPS_NO_PROPOGATED_SPLITS.mgf.tar.gz; then
    t2=$(stat -c %y $LIBRARY_FILE)
    if [ "$t1" != "$t2" ]; then 
        echo "Download Successful! A new version of the ALL_GNPS_NO_PROPOGATED library was retrieved! "
        # if wget was perfomed, untar the downloaded file
        tar -xvf ALL_GNPS_NO_PROPOGATED_SPLITS.mgf.tar.gz
        # join the split content
        cat ALL_GNPS_NO_PROPOGATED_*.mgf > ALL_GNPS_NO_PROPOGATED.mgf
        # remove the split parts
        rm ALL_GNPS_NO_PROPOGATED_*.mgf  
        # lets create a summary for the library file - in this case only the ALL_GNPS_NO_PROPOGATED.mgf is used, parse the ions headers into a table
        echo "Creating a summary of the ALL_GNPS_NO_PROPOGATED library!"
        cd ../..
        #python GNPS2_LibrarySearch_Workflow/library_summary.py libraries/GNPS-LIBRARY.mgf data/library_summary_GNPS-LIBRARY.tsv
        python GNPS2_LibrarySearch_Workflow/library_summary.py data/libraries/ALL_GNPS_NO_PROPOGATED.mgf data/library_summary_ALL_GNPS_NO_PROPOGATED.tsv
        echo "Enriching the ALL_GNPS_NO_PROPOGATED library summary with other annotations!"
        # execute a join between the library summary and the library summary enriched with annotations (hardcoded)
        #python GNPS2_LibrarySearch_Workflow/library_summary_merge_annotations.py data/library_summary_GNPS-LIBRARY.tsv data/library_summary_GNPS-LIBRARY_annotations.tsv data/library_summary_GNPS-LIBRARY_enriched.tsv
        python GNPS2_LibrarySearch_Workflow/library_summary_merge_annotations.py data/library_summary_ALL_GNPS_NO_PROPOGATED.tsv data/library_summary_ALL_GNPS_NO_PROPOGATED_annotations.tsv data/library_summary_ALL_GNPS_NO_PROPOGATED_enriched.tsv      
    else
        echo "The ALL_GNPS_NO_PROPOGATED library is already present in its latest version."; fi
    fi   
else
    echo "Download Failed for ALL_GNPS_NO_PROPOGATED library =("
    echo "Check if there is enough space in your disk, at least 7Gb are necessary for downloading and joining the splitted parts. At the end, 3Gb will be used to store the full library."
fi

#wget https://external.gnps2.org/gnpslibrary/GNPS-SELLECKCHEM-FDA-PART1.mgf
#cd ../ && mkdir spectra && cd spectra
#wget --output-document=isa_9.mzML "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000084030/ccms_peak/isa_9.mzML&forceDownload=true"
