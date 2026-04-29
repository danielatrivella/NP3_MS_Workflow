# @Crisfbazz
# adapted from getGNPS_library_annotations.py to enrich a library summary table with online data
# this is intend to be executed once by the adm with internet connection (online)
# then the library search result can be enriched offline with all these retrieved info - workaround to execute np3 offline and remove dependency to the GNPS server for annotation enrichment
# this should be updated by the adm everytime the Library updates (new spectrumIDs are added) - a warning will show up in the joining step from NP3 setup
python ../GNPS2_LibrarySearch_Workflow/getGNPS_library_summary_annotations.py library_summary_ALL_GNPS_NO_PROPOGATED.tsv library_summary_ALL_GNPS_NO_PROPOGATED_annotations.tsv --numthreads 300
# compress resulting table
tar -cf - library_summary_ALL_GNPS_NO_PROPOGATED_annotations.tsv | gzip -9 > library_summary_ALL_GNPS_NO_PROPOGATED_annotations.tsv.tar.gz
