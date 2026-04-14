# @Crisfbazz
# Adapted from https://github.com/Wang-Bioinformatics-Lab/LibrarySearch_Workflow/blob/master/data/get_data.sh
# commit 4284470
# Following the flow present in the GNPS2 LibrarySearch_Workflow
# adapted here to be executed offline and using only the GNPS-library.mgf library file

# mkdir if not created yet
mkdir -p libraries && cd libraries
# download the library file if a new version is present in the server
wget -N https://external.gnps2.org/gnpslibrary/GNPS-LIBRARY.mgf

#wget https://external.gnps2.org/gnpslibrary/GNPS-SELLECKCHEM-FDA-PART1.mgf
#cd ../ && mkdir spectra && cd spectra
#wget --output-document=isa_9.mzML "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000084030/ccms_peak/isa_9.mzML&forceDownload=true"

# lets create a summary for the library file - in this case only the GNPS-library.mgf will be downloaded and used
cd ../
python ../GNPS2_LibrarySearch_Workflow/library_summary.py libraries/GNPS-LIBRARY.mgf library_summary.tsv