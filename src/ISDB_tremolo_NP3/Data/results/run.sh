# Remember, never use spaces or special characters in your files
#
# By default:
# - Uncompress tremolo and put it in a tremolo folder
# - Put the mgf files of the database in a dbs folder
# - Copy the UNPD_DB.csv file in the dbs folder
# - Copy treat.py in the dbs folder
# - Put this script in a "results" folder (this doesn't change much)
# - Go to the results folder and run: bash ../dbs/run.sh your_clustered_spectra_file.mgf your_cytoscape_attribute_file.out your_output.out


# Give the path of the libraries relative to where you run this script
# They should be separated by spaces
EXISTING_LIBRARY_MGF="../dbs/UNPD_ISDB_R_p01.mgf ../dbs/UNPD_ISDB_R_p02.mgf ../dbs/UNPD_ISDB_R_p03.mgf ../dbs/UNPD_ISDB_R_p04.mgf ../dbs/UNPD_ISDB_R_p05.mgf ../dbs/UNPD_ISDB_R_p06.mgf ../dbs/UNPD_ISDB_R_p07.mgf ../dbs/UNPD_ISDB_R_p08.mgf ../dbs/UNPD_ISDB_R_p09.mgf"
# Give the path to the CSV file containing the description of the database
DB_DESC="../dbs/UNPD_DB.csv"

# Set the tolerance to be used
TOLERANCE=0.025
# Score threshold
SCORE_THRESHOLD=0.2
# Top K results
TOP_K_RESULTS=3


if [ $# -lt 3 ]
then
    echo "Usage: $0 <spectrafile as mgf> <cytoscape_annotation> <output_file_annotation>"
    exit
fi

if [ ! -f $1 ]
then
   echo "The spectra file $1 doesn't exist"
   exit
fi

if [ ! -f $2 ]
then
    echo "The cytoscape annotation file $2 doesn't exist"
    exit
fi


if [ "$2" == "$3" ]
then
    echo "The cytoscape input and output cannot be the same file"
    exit
fi
# Name of the output file that you will use in your graph software
OUTPUT_FILE=$3


# Treating script
TREAT="../dbs/treat.py"
# Give the tremolo path
TREMOLO_PATH="../tremolo/"
# Name of the tremolo output
RESULT_FILE="./Results_tremolo.out"
# Name of the actual annotation file
INPUT_FILE=$2

# Merging database

if [ ! -f ../dbs/UNPD_DB.csv ]
then
    echo "Merging database"
    OLDPATH=`pwd`
    cd ../dbs
    ./merge_db.sh
    cd $OLDPATH
fi

# Convert the spectra file as a pkllib file
echo "Converting the spectra file"
SPECT_FILE=`basename $1`.pklbin
if [ -f $SPECT_FILE ]
then
    rm $SPECT_FILE
fi

$TREMOLO_PATH/convert $1 $SPECT_FILE

if [ ! -f $SPECT_FILE ]
then
    echo " Conversion failed"
    exit
fi



# Create the parameters file
echo "Creating the parameters file"
if [ -f scripted.params ]
then
    rm scripted.params    
fi

cat > scripted.params <<EOF 
EXISTING_LIBRARY_MGF=$EXISTING_LIBRARY_MGF

searchspectra=./$SPECT_FILE

RESULTS_DIR=$RESULT_FILE

tolerance.PM_tolerance=$TOLERANCE

search_decoy=0

SCORE_THRESHOLD=$SCORE_THRESHOLD
TOP_K_RESULTS=$TOP_K_RESULTS

NODEIDX=0
NODECOUNT=1

SEARCHABUNDANCE=0
SLGFLOADMODE=1
EOF

echo "Running tremolo (the output is logged in tremolo.log)"
$TREMOLO_PATH/main_execmodule ExecSpectralLibrarySearch ./scripted.params &> tremolo.log

#echo "Running the treatment script"
#bash run.sh ~/Documents/CNPEM/ms_entity_annotation/Data_Collections/Bra346/Bra_RT1_2_test/outs/#Bra_RT1_2_test/mgf/Bra_RT1_2_test_all_clean.mgf 
#run.sh results.out

python2.7 $TREAT $RESULT_FILE $INPUT_FILE $DB_DESC $OUTPUT_FILE
#python2.7 ../dbs/treat.py ./Results_tremolo.out #Bra_RT1_2_test_peak_area_count_clean_annotated.csv ../dbs/UNPD_DB.csv bra_result.out

#echo "You can now use the file $OUTPUT_FILE as annotation in your graph program"
