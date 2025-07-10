# script to join the ivamns from different jobs being joined
# first the original msclusterIDs of each IVAMN is mapped to their msclusterID + _ + JOB_CODE
# then the original msclusterIDs with job code will be replaced by the new msclusterID of the joined job
# using the joinedJobsIDs of the clean table as reference
# duplicated edges are merged at the end and
# finally the protonated script is executed for the final IVAMN and the final list of protonated_representatives
# is merged with the original protonated that have in-degree > 0 in the final net.
# the joined annotations will be written to the clean table as columns in another script
# If noise cutoff is informed, also remove original spectra below the minimum basePeakInt from the final clean table

import os, sys
import pandas as pd
import numpy as np
from mn_annotations_assign_protonated_representative import mn_annotation_find_protonated

#  Root Mean Square Error (RMSE) is a metric used to evaluate the difference between predicted and actual values
# used here to compute the rt error between annotated nodes
def rmse(predicted, actual):
    return np.sqrt(np.mean((predicted - actual)**2))

# remove edges with the invalid ann below the sim cutoff
# check if the selected edges to remove also contain another type of annotation, if yes only remove the fragment type
#not_valid_ann = neutral losses < 0.2, fragments < 0.2 or isotopes < 0.75 with low cosine
#not_valid_ann="(\\[M\\+H-NH3]\\+)|(\\[M\\+H-H2O]\\+)|(\\[M\\+H-NH3-H2O]\\+)" or "(fragments)" or "(\\[M\\+1\\]\\+)"
def remove_not_valid_ann_from_ivamn(ivamn, not_valid_ann, cosine_cutoff=0.0):
    # remove edges with only the invalid annotation and below the cutoff
    ivamn = ivamn.loc[~(ivamn.annotation.str.fullmatch(not_valid_ann) & (ivamn.cosine < cosine_cutoff)), :].copy()
    # remove the invalid annotation from edges below the cutoff and with more than one annotation
    # check if any of the groups is present, them extract them and remove them from the annotations of each edge
    ann_with_not_valid = (ivamn.annotation.str.contains(not_valid_ann.replace("(","").replace(")",""), regex=True) &
                          (ivamn.cosine < cosine_cutoff))
    if ann_with_not_valid.any():
        # extract the matched string
        matched_not_valid_ann = ivamn.annotation[ann_with_not_valid].str.extract(not_valid_ann)[0].values
        ivamn.loc[ann_with_not_valid, "annotation"] = np.asarray([ann.replace(matched_not_valid_ann[i],"").lstrip(";").rstrip(";").replace(";;",";")
                                                       for i, ann in enumerate(ivamn.annotation[ann_with_not_valid].values)])
    return ivamn


# read the join original jobs metadata
# for each original job being joined, read its IVAMN and map its msclusterID to the joined msclusterID
# then concatenate the IVAMNS and merge duplicated edges
def join_jobs_ivamns(output_path, max_chunk=3000, noise_cutoff=False):
    if not os.path.isdir(output_path):
        sys.exit("ERROR. The provided output path directory '"+output_path+"' does not exists.")
    output_name = os.path.basename(output_path)

    print("* Joining the Original IVAMNs of the Joined Jobs * ")
    # read the clean quantifications of the output_name job
    clean_counts_path = os.path.join(output_path, "count_tables", "clean", output_name+"_peak_area_clean.csv")
    if not os.path.isfile(clean_counts_path):
        sys.exit("ERROR. The clean quantification of the current joining jobs in terms of peak area '" +
                 clean_counts_path + "' does not exists.")
    clean_counts = pd.read_csv(clean_counts_path, low_memory=False, usecols=['msclusterID', 'joinedJobsIDs',
                                                                             'mzConsensus', 'rtMean', 'rtMin', 'rtMax',
                                                                             'basePeakInt'])
    clean_counts.sort_values(by=['msclusterID'], inplace=True, ignore_index=True)
    # compute a noise cutoff as the minimum basePeakInt after union if noise cutoff is enabled
    # use the min basePeakInt in the final clean table as the basePeakInt cutoff
    if noise_cutoff:
        basePeakInt_noise_cutoff = clean_counts.basePeakInt.min()
        print("  - Noise cut-off was enabled and set to ", str(basePeakInt_noise_cutoff))
    else:  # noise cutoff is disabled
        basePeakInt_noise_cutoff = 0

    # read the default join original jobs metadata, from it extract the Job name, code and path
    original_jobs_metadata_path = os.path.join(output_path, "../..", "original_jobs_METADATA_JOIN.csv")
    if not os.path.isfile(original_jobs_metadata_path):
        sys.exit("ERROR. The metadata table with the original joined jobs '" + original_jobs_metadata_path +
                 "' does not exists.")
    jobs_metadata = pd.read_csv(original_jobs_metadata_path)

    # check if the similarity table exists
    sim_table_path = os.path.join(output_path, "molecular_networking", "similarity_tables",
                                  "similarity_table_" + output_name + "_clean.csv")
    if not os.path.isfile(sim_table_path):
        sys.exit("ERROR. The clean similarity table of the current job being joined with name '" + output_name +
                 "' does not exists." +
                 " The following path constructed with the provided output_path is not valid: '" + sim_table_path +
                 "'.")

    if max_chunk < 2:
        print("WARNING: the max_chunk parameter have a very low value smaller than 2. ",
              "It will be set to 3000 as the default value.")
        max_chunk = 3000

    # for each job code plus one iteration for integration of current output_name (final one):
    # read its ivamn, map the original msclusterID of each job to the joined jobs IDs
    # use the ivamns attributes tables first? then map to the ivamn net?
    # from the attributes tables retrieve the multicharge_ion and isotope_ion info and protonated_representative (old)
    # recompute the protonated again and compare
    # retrieve similarity from clean sim table
    for i in range(jobs_metadata.shape[0]+1):
        if i == jobs_metadata.shape[0]:
            # final integration step, retrieve the temporary ivamn att table from the output path
            # read the ivamn attribute table
            job_ivamn_att_path = os.path.join(output_path, "molecular_networking",
                                              output_name + "_ivamn_attributes_tmp.csv")
            job_code = output_name
            print("    - Integrating and reducing the mapped IVAMNs for current joined job " + job_code)
        else:
            job_code = jobs_metadata.JOB_CODE[i]
            # read the ivamn attribute table
            job_ivamn_att_path = os.path.join(jobs_metadata.JOB_PATH[i], "molecular_networking",
                                                 jobs_metadata.JOBNAME[i]+"_ivamn_attributes.csv")
            print("    - Mapping and reducing the IVAMN of job code " + job_code)
        if not os.path.isfile(job_ivamn_att_path):
            sys.exit("ERROR. The IVAMN attribute table of the original job with code '"+job_code +
                     "' does not exists. The following path extracted from the joined jobs metadata is not valid: '" +
                     job_ivamn_att_path + "'.")
        job_ivamn_att = pd.read_csv(job_ivamn_att_path, low_memory=False,
                                    usecols=['msclusterID', 'multicharge_ion', 'isotope_ion',
                                             'protonated_representative'])
        # also read basePeakInt from the clean table of the original file and remove spectra by the noise cutoff if enabled
        # remove nodes and edges between removed nodes using the noise cutoff
        if i < jobs_metadata.shape[0]:
            # if this is not the final integration step
            # remove msclusterIDs with basePeakInt below the cutoff if enabled
            # and then map the msclusterID_job to the new joinedJobs msclusterID

            if noise_cutoff:
                job_clean_path = os.path.join(jobs_metadata.JOB_PATH[i], "count_tables", "clean",
                                                  jobs_metadata.JOBNAME[i] + "_peak_area_clean_ann.csv")
                if not os.path.isfile(job_clean_path):
                    sys.exit("ERROR. The clean table of the original job with code '" + job_code +
                             "' does not exists. The following path extracted from the joined jobs metadata is not valid: '" +
                             job_clean_path + "'.")
                job_clean = pd.read_csv(job_clean_path, low_memory=False, usecols=['msclusterID', 'basePeakInt'])
                valid_msclusterIDs = job_clean.msclusterID[job_clean.basePeakInt >= basePeakInt_noise_cutoff]
                del job_clean
                # filter job_ivamn_att and then filter ivamn net edges
                job_ivamn_att = job_ivamn_att.loc[job_ivamn_att.msclusterID.isin(valid_msclusterIDs),]

            job_ivamn_att["msclusterID_job"] = job_ivamn_att.msclusterID.astype(str)+ "_" + job_code
            try:
                job_ivamn_att["msclusterID_new"] = [clean_counts.msclusterID.values[
                     clean_counts.joinedJobsIDs.str.contains("^"+id + "$|^" + id + ";|;" + id + "$|;" + id + ";", regex=True)][0]
                 for id in job_ivamn_att.msclusterID_job.values]
            except IndexError:
                # error if any id is not found and indexing the first position of values return an error
                missing_ids = np.asarray([clean_counts.joinedJobsIDs.str.contains("^" + id + "$|^" + id + ";|;" + id + "$|;" + id + ";", regex=True).any()
                 for id in job_ivamn_att.msclusterID_job.values])
                sys.exit("ERROR. Index error when mapping the original job '" + job_code +
                         "' with the joined jobs IDs. The following original IDs are missing from the joined clean table: " +
                         ','.join(job_ivamn_att.msclusterID_job[np.where(~missing_ids)[0]].values) + ". Wrong mapping.")
            # set the original msclusterID as index of the att table
            job_ivamn_att = job_ivamn_att.set_index('msclusterID')
        #
        if i == jobs_metadata.shape[0]:
            # final integration step, retrieve the temporary ivamn net from the output path
            # read the ivamn and fix the msclusterID columns types to int
            job_ivamn_path = os.path.join(output_path, "molecular_networking",
                                          output_name + "_ivamn_tmp.selfloop")
        else:
            # read the ivamn and fix the msclusterID columns types to int
            job_ivamn_path = os.path.join(jobs_metadata.JOB_PATH[i], "molecular_networking",
                                          jobs_metadata.JOBNAME[i] + "_ivamn.selfloop")
        if not os.path.isfile(job_ivamn_path):
            sys.exit("ERROR. The IVAMN of the original job with code '" + job_code + "' does not exists." +
                     " The following path extracted from the joined jobs metadata is not valid: '" + job_ivamn_path +
                     "'.")
        job_ivamn = pd.read_csv(job_ivamn_path,
                                dtype={'msclusterID_source': np.int64, 'msclusterID_target': np.int64},
                                low_memory=False, usecols= ['msclusterID_source', 'msclusterID_target', 'annotation',
                                                            'mzError', 'rtError'])
        # remove self loops, only map the not isolated nodes - selfloops will be added at the end
        job_ivamn = job_ivamn.loc[(job_ivamn.msclusterID_source != job_ivamn.msclusterID_target), :]
        #
        if i < jobs_metadata.shape[0]:
            # if this is not the final integration step
            # map the msclusterID_job to the new joinedJobs msclusterID
            # use the ivamn attribute table with the original msclusterIDs as index to map the source and target IDs to
            # # the msclusterID_new in the joined jobs
            # first remove the edges between not valid nodes
            if noise_cutoff:
                # filter ivamn net edges between valid msclusterIDs
                job_ivamn = job_ivamn.loc[(job_ivamn.msclusterID_source.isin(valid_msclusterIDs) &
                                           job_ivamn.msclusterID_target.isin(valid_msclusterIDs)), :]
            if job_ivamn.shape[0] > 0:
                # if any valid edges left, proceed for the mapping
                job_ivamn["msclusterID_source_new"] = job_ivamn_att.loc[
                    job_ivamn.msclusterID_source.values, "msclusterID_new"].values
                job_ivamn["msclusterID_target_new"] = job_ivamn_att.loc[
                    job_ivamn.msclusterID_target.values, "msclusterID_new"].values
                job_ivamn = job_ivamn.drop(["msclusterID_source", "msclusterID_target"], 1)
                # rename the new msclusterID columns
                job_ivamn.rename(columns={'msclusterID_source_new': 'msclusterID_source',
                                          'msclusterID_target_new': 'msclusterID_target'}, inplace=True)
        duplicated_joined_edges = np.where(job_ivamn.loc[:,["msclusterID_source","msclusterID_target"]].duplicated(keep=False))[0]
        first_duplicated_edges = duplicated_joined_edges[
            np.where(~job_ivamn.loc[duplicated_joined_edges, ["msclusterID_source","msclusterID_target"]].duplicated(keep="first"))[0]]
        # for each duplicated edges, concatenate the unique annotations and compute the mean mzError
        # the mzError was computed based on the type of the annotation, so here the mean is used to ease the computations
        # the rtError will be recomputed at the end
        for dup_edge in first_duplicated_edges:
            dup_source_id, dup_target_id = job_ivamn.loc[dup_edge, ["msclusterID_source", "msclusterID_target"]]
            select_duplicated_rows = ((job_ivamn.loc[duplicated_joined_edges, "msclusterID_source"] == dup_source_id) &
                (job_ivamn.loc[duplicated_joined_edges, "msclusterID_target"] == dup_target_id))
            unique_annotations_concat = ';'.join(job_ivamn.loc[duplicated_joined_edges,:].loc[select_duplicated_rows,
                "annotation"].unique())
            mzError_mean = job_ivamn.loc[duplicated_joined_edges, :].loc[select_duplicated_rows, "mzError"].mean()
            job_ivamn.loc[dup_edge, "annotation"] = unique_annotations_concat
            job_ivamn.loc[dup_edge, "mzError"] = np.round(mzError_mean, 3)
        #
        # remove the duplicated edges, except the first duplicated ones
        job_ivamn = job_ivamn.drop(np.setdiff1d(duplicated_joined_edges, first_duplicated_edges), 0)
        # order the columns and store the ivamn temporary in the output_path
        job_ivamn = job_ivamn[['msclusterID_source','msclusterID_target', 'annotation', 'mzError', 'rtError']]
        if i == 0:
            job_ivamn.to_csv(os.path.join(output_path, "molecular_networking", output_name+"_ivamn_tmp.selfloop"),
                             index=False)
        else:
            if i < jobs_metadata.shape[0]:
                # append to temporary ivamn
                job_ivamn.to_csv(os.path.join(output_path, "molecular_networking", output_name + "_ivamn_tmp.selfloop"),
                                 index=False, mode='a', header=False)
            else:
                # if this is the final integration step, overwrite the tmp file
                job_ivamn.sort_values(by=["msclusterID_source", "msclusterID_target"], inplace=True, ignore_index=True)
                job_ivamn.to_csv(os.path.join(output_path, "molecular_networking", output_name + "_ivamn_tmp.selfloop"),
                                 index=False)
        # rm current ivamn
        del job_ivamn
        #
        if i < jobs_metadata.shape[0]:
            # if this is not the final integration step, rename the new msclusterID mapped column
            job_ivamn_att.rename(columns={'msclusterID_new': 'msclusterID'}, inplace=True)
        # remove duplicates from the att table and store it temporary in the output path
        job_ivamn_att = job_ivamn_att.reset_index(drop=True)
        duplicated_joined_nodes = np.where(job_ivamn_att.loc[:, "msclusterID"].duplicated(keep=False))[0]
        first_duplicated_nodes = duplicated_joined_nodes[np.where(
            ~job_ivamn_att.loc[duplicated_joined_nodes, "msclusterID"].duplicated(keep="first"))[0]]
        # for each duplicated node, aggregate the multicharge_ion, the isotope_ion and the protonated_representative columns
        # make them receive the max values among the joined ids
        for dup_node in first_duplicated_nodes:
            dup_id = job_ivamn_att.loc[dup_node, "msclusterID"]
            select_duplicated_rows = (job_ivamn_att.loc[duplicated_joined_nodes, "msclusterID"] == dup_id)
            job_ivamn_att.loc[dup_node, ["multicharge_ion", "isotope_ion", "protonated_representative"]] = job_ivamn_att.loc[duplicated_joined_nodes,:].loc[
                select_duplicated_rows, ["multicharge_ion", "isotope_ion", "protonated_representative"]].max()

        #
        # remove the duplicated nodes, except the first duplicated ones
        job_ivamn_att = job_ivamn_att.drop(np.setdiff1d(duplicated_joined_nodes, first_duplicated_nodes), 0)
        # order columns and store the att table
        job_ivamn_att = job_ivamn_att[['msclusterID', 'multicharge_ion', 'isotope_ion', 'protonated_representative']]
        # check if all the joined msclusterIDs are in the clean counts table ids, trought error if not
        if not job_ivamn_att.msclusterID.isin(clean_counts.msclusterID).all():
            sys.exit("ERROR. The joined IVAMN attribute table of the original job with code '" + job_code +
                     "' was not correctly created. There is a final msclusterID that do not appear in the joined clean " +
                     " quantification table (problematic IDs: " +
                     ','.join(job_ivamn_att.msclusterID[~job_ivamn_att.msclusterID.isin(clean_counts.msclusterID)]) +
                     "). Something went wrong when joining the IVAMN att table and reducing redundancies.")
        if i == 0:
            job_ivamn_att.to_csv(os.path.join(output_path, "molecular_networking", output_name + "_ivamn_attributes_tmp.csv"),
                             index=False)
        else:
            if i < jobs_metadata.shape[0]:
                # append to temporary ivamn
                job_ivamn_att.to_csv(os.path.join(output_path, "molecular_networking", output_name + "_ivamn_attributes_tmp.csv"),
                                 index=False, mode='a', header=False)
            else:
                # if this is the final integration step, check if the number of rows match the nrows of the clean counts
                # # and overwrite the tmp file
                if job_ivamn_att.shape[0] != clean_counts.shape[0]:
                    sys.exit("ERROR. The joined IVAMN attribute table of the joined job with name '" + job_code +
                             "' was not correctly created. The final number of rows of its joined IVAMN attribute table ("+
                             str(job_ivamn_att.shape[0])+") " +
                             "does not match with the number of rows in its clean quantification table (" +
                             str(clean_counts.shape[0]) +
                             "). Something went wrong when integrating the joined IVAMNs att table and reducing redundancies.")
                # order the att table by the msclusterID column
                job_ivamn_att.sort_values(by=['msclusterID'], inplace=True, ignore_index=True)
                # order columns and store the att table
                job_ivamn_att = job_ivamn_att[['msclusterID', 'multicharge_ion', 'isotope_ion',
                                               'protonated_representative']]
                job_ivamn_att.rename(columns={'protonated_representative': 'protonated_representative_old'}, inplace=True)
                job_ivamn_att.to_csv(job_ivamn_att_path, index=False)
        # rm current ivamn att table
        del job_ivamn_att

    # At the end, after the joined ivamn integration
    # recompute rtError and retrieve the new cosine for each available edge
    # also recompute number of common samples
    # Finally, call the protonated representative selection, which also retrieves the componentIndex
    job_ivamn = pd.read_csv(job_ivamn_path,
                            dtype={'msclusterID_source': np.int64, 'msclusterID_target': np.int64},
                            low_memory=False, usecols=['msclusterID_source', 'msclusterID_target', 'annotation',
                                                       'mzError', 'rtError'])
    # sort the job_ivamn by the minimum msclusterID - create column msclusterID_min - then use this to retrieve
    # the cosine from the pairwise table by chunks,
    # when the current min > max_chunk, read the next chunk of the pairwise table
    # always read the cosine from the min msclusterID x max msclusterID in the edge
    job_ivamn[['msclusterID_min', 'msclusterID_max']] = job_ivamn[['msclusterID_source', 'msclusterID_target']].apply(
        lambda x: [min(x), max(x)], axis=1, result_type="expand")
    job_ivamn.sort_values(by=['msclusterID_min'], inplace=True, ignore_index=True)
    # read the clean_counts again to retrieve the quantification columns by samples - used to compute the number of common samples
    clean_counts = pd.read_csv(clean_counts_path, low_memory=False)
    # filter columns 'msclusterID', 'joinedJobsIDs', 'mzConsensus', 'rtMean', 'rtMin', 'rtMax' and the ones that ends with _area
    clean_counts = clean_counts[np.concatenate((np.asarray(['msclusterID', 'joinedJobsIDs', 'mzConsensus',
                                                            'rtMean', 'rtMin', 'rtMax']),
                                clean_counts.columns[clean_counts.columns.str.endswith("_area")].values))]
    clean_counts.sort_values(by=['msclusterID'], inplace=True, ignore_index=True)
    # get the count_columns indexes to perform and xand operation for each edge
    counts_cols_idx = np.where(clean_counts.columns.str.endswith("_area"))[0]
    # read the similarity table by chunk - sim_table_path e max_chunk
    sim_table = pd.read_csv(sim_table_path, low_memory=False, nrows=max_chunk)
    current_chunk = 0  # set the first line of the sim_table
    # for each edge in the new joined ivamn, retrieve and recompute the rtError, cosine and numCommonSamples
    job_ivamn['cosine'] = 0
    job_ivamn['numCommonSamples'] = 0
    print("    - Recomputing the edges' attributes of the joined IVAMN ")
    for i in range(job_ivamn.shape[0]):
        # the idx position of the connected nodes in the clean_counts and similarity table - same ordering
        source_i = np.where(clean_counts.msclusterID == job_ivamn.msclusterID_min[i])[0]
        target_i = np.where(clean_counts.msclusterID == job_ivamn.msclusterID_max[i])[0]
        # recompute the rtError
        rtMean_i, rtMin_i, rtMax_i = clean_counts.loc[source_i, ['rtMean', 'rtMin', 'rtMax']].values[0]
        rtMean_i_j, rtMin_i_j, rtMax_i_j = clean_counts.loc[target_i, ['rtMean', 'rtMin', 'rtMax']].values[0]
        # compute the rtError for blank and not blank nodes
        if (((rtMin_i == 0) | (rtMin_i_j == 0)) & ((rtMax_i == 1000000) | (rtMax_i_j == 1000000))):
            # this is an edge between a blank node, use only the rtMean to compute the rtError
            rtError_i = np.round(rmse(rtMean_i, rtMean_i_j), 2)
        else:
            # this is an edge between not blank nodes, use the rtMin and rtMax to compute the rtError
            rtError_i = np.round(rmse(np.asarray([rtMin_i,rtMax_i]), np.asarray([rtMin_i_j,rtMax_i_j])), 2)
        # set the new rtError to the ivamn net
        job_ivamn.loc[i, "rtError"] = rtError_i
        # recompute the numCommonSamples
        num_common_samples = (clean_counts.iloc[np.concatenate((source_i,target_i)),
                                                counts_cols_idx] > 0).all(0).sum()
        job_ivamn.loc[i, "numCommonSamples"] = np.int64(num_common_samples)
        # retrieve the new cosine
        if source_i >= current_chunk+max_chunk:
            # read next chunk of the sim_table
            sim_table = pd.read_csv(sim_table_path, low_memory=False, nrows=max_chunk, skiprows=source_i[0]+1,
                                    header=None, names=sim_table.columns)
            current_chunk = source_i[0]
        pairwise_sim = sim_table.iloc[source_i-current_chunk, target_i + 1].values[0]
        job_ivamn.loc[i, "cosine"] = np.round(pairwise_sim, 3)
        # round mzError again
        job_ivamn.loc[i, "mzError"] = np.round(job_ivamn.loc[i, "mzError"], 3)

    print("    - Removing fragment's and neutral losses annotations with similarity < 0.2 and isotope annotations with similarity < 0.75")
    # remove annotations of fragments with similarity values < 0.2
    # remove annotations of isotopes with similarity values < 0.75, of neutral losses with similarity values < 0.2
    # defined the regex to select the not valid annotations using groups: ()
    job_ivamn = remove_not_valid_ann_from_ivamn(job_ivamn, not_valid_ann="(fragment)", cosine_cutoff=0.2)
    job_ivamn = remove_not_valid_ann_from_ivamn(job_ivamn, not_valid_ann="(\\[M\\+H-NH3]\\+)|(\\[M\\+H-H2O]\\+)|(\\[M\\+H-NH3-H2O]\\+)",
                                                cosine_cutoff=0.2)
    job_ivamn = remove_not_valid_ann_from_ivamn(job_ivamn,
                                                not_valid_ann="(\\[M\\+1\\]\\+)",
                                                cosine_cutoff=0.75)
    # order columns of the final ivamn
    job_ivamn = job_ivamn[['msclusterID_source', 'msclusterID_target', 'cosine', 'annotation', 'mzError',
                           'rtError', 'numCommonSamples']]
    # sort by msclusterID_source and msclusterID_target and save final IVAMN
    job_ivamn.sort_values(by=["msclusterID_source", "msclusterID_target"], inplace=True, ignore_index=True)
    print("    - Adding selfloop edges to the IVAMN")
    # add to the IVAMN the selfloop edges connecting the isolated nodes
    isolated_nodes = clean_counts.msclusterID[~(clean_counts.msclusterID.isin(job_ivamn.msclusterID_source) |
                                                clean_counts.msclusterID.isin(job_ivamn.msclusterID_target))].values
    job_ivamn = pd.concat([job_ivamn, pd.DataFrame({'msclusterID_source': isolated_nodes,
                                                   'msclusterID_target': isolated_nodes, 'cosine': 1.00,
                                                   'annotation': '', 'mzError': 0.0, 'rtError': 0.0,
                                                   'numCommonSamples': ''})])
    # save the final IVAMN
    job_ivamn_path = os.path.join(output_path, "molecular_networking", output_name + "_ivamn.selfloop")
    job_ivamn.to_csv(job_ivamn_path, index=False)
    print("    - Adding the multicharge_ion and isotope_ion to the clean count tables")
    # join the multicharge_ion and isotope_ion cols from tmp ivamn att to the clean table
    # store in the clean table _ann
    job_ivamn_att = pd.read_csv(job_ivamn_att_path, low_memory=False,
                                usecols=['msclusterID', 'multicharge_ion', 'isotope_ion', 'protonated_representative_old'])
    job_ivamn_att.sort_values(by=['msclusterID'], inplace=True, ignore_index=True)
    clean_counts = pd.read_csv(clean_counts_path, low_memory=False)
    clean_counts.sort_values(by=['msclusterID'], inplace=True, ignore_index=True)
    clean_counts[['multicharge_ion', 'isotope_ion']] = job_ivamn_att[['multicharge_ion', 'isotope_ion']]
    clean_counts.to_csv(os.path.join(output_path, "count_tables", "clean", output_name + "_peak_area_clean_ann.csv"),
                        index=False)
    # also merge in the spectra count
    clean_counts_path = os.path.join(output_path, "count_tables", "clean", output_name + "_spectra_clean.csv")
    if not os.path.isfile(clean_counts_path):
        sys.exit("ERROR. The clean quantification of the current joining jobs in terms of number of spectra '" +
                 clean_counts_path + "' does not exists.")
    clean_counts = pd.read_csv(clean_counts_path, low_memory=False)
    clean_counts.sort_values(by=['msclusterID'], inplace=True, ignore_index=True)
    clean_counts[['multicharge_ion', 'isotope_ion']] = job_ivamn_att[['multicharge_ion', 'isotope_ion']]
    clean_counts.to_csv(os.path.join(output_path, "count_tables", "clean", output_name + "_spectra_clean_ann.csv"),
                        index=False)
    # set the path to the peak area count again
    clean_counts_path = os.path.join(output_path, "count_tables", "clean", output_name + "_peak_area_clean_ann.csv")
    print("\n* Assigning the putative [M+H]+ spectra representatives in the joined IVAMN and merging with original [M+H]+ that have in-degree > 0 *")
    # call find protonated - this will create the final ivamn att table and add componentIndex to the IVAMN
    mn_annotation_find_protonated(job_ivamn_path, clean_counts_path)
    # merge the final ivamn att table
    # merge the old protonated with the new list of protonated_representative
    # add as protonated the old protonated nodes with in-degree > 0 in the final ivamn
    # this removes old protonated that were selfloops in the original jobs and now are in a cluster with at least
    # # one in_degree, otherwise it is not be added
    job_ivamn_att_final = pd.read_csv(os.path.join(output_path, "molecular_networking", output_name + "_ivamn_attributes.csv"),
                                      low_memory=False)
    job_ivamn_att_final.sort_values(by=['msclusterID'], inplace=True, ignore_index=True)
    job_ivamn_att_final.loc[((job_ivamn_att.protonated_representative_old == 1) & (job_ivamn_att_final.in_degree > 0)),
                            "protonated_representative"] = 1
    print("Number of protonated_representative nodes (with isolated nodes) after merge with old protonated",
          "that have in-degree > 0:", str(sum(job_ivamn_att_final.protonated_representative == 1)))
    job_ivamn_att_final.to_csv(os.path.join(output_path, "molecular_networking", output_name + "_ivamn_attributes.csv"),
                               index=False)
    del job_ivamn_att
    # update the protonated_representative in the clean tables
    clean_counts = pd.read_csv(clean_counts_path, low_memory=False)
    clean_counts.sort_values(by=['msclusterID'], inplace=True, ignore_index=True)
    clean_counts['protonated_representative'] = job_ivamn_att_final['protonated_representative']
    clean_counts.to_csv(clean_counts_path, index=False)
    clean_counts_path = os.path.join(output_path, "count_tables", "clean", output_name + "_spectra_clean_ann.csv")
    clean_counts = pd.read_csv(clean_counts_path, low_memory=False)
    clean_counts.sort_values(by=['msclusterID'], inplace=True, ignore_index=True)
    clean_counts['protonated_representative'] = job_ivamn_att_final['protonated_representative']
    clean_counts.to_csv(clean_counts_path, index=False)
    del clean_counts
    # remove the tmp ivamn att table and network files
    os.remove(os.path.join(output_path, "molecular_networking",
                                          output_name + "_ivamn_tmp.selfloop"))
    os.remove(os.path.join(output_path, "molecular_networking",
                                              output_name + "_ivamn_attributes_tmp.csv"))


if __name__ == "__main__":
    max_chunk = 3000
    noise_cutoff = False
    if len(sys.argv) > 1:
        # print(sys.argv)
        output_path = sys.argv[1]
        if len(sys.argv) > 2:
            max_chunk = int(sys.argv[2])
            if len(sys.argv) > 3:
                noise_cutoff = bool(float(sys.argv[3]))
    else:
        print("Error: One argument must be supplied to join the original IVAMNs of jobs being joined inside the join_jobs command flow:\n",
        " 1 - output_path: Path to the final output directory of the job currently being joined, ",
        "named with the output_name inside the 'outs' directory;\n",
        " 2 - max_chunk: Maximum number of rows of the pairwise similarity table to be loaded and ",
        "processed in a chunk at the same time. In case of memory issues this value should be decreased (default: 3000);\n",
        " 3 - noise_cutoff: A boolean True or False indicating if the noise cutoff was enable in cleaning, if yes ",
        "remove spectra from original samples with a basePeakInt < minimum basePeakInt in the final joined clean table; ",
        "otherwise do nothing (default: False).\n")
        sys.exit(1)
    # call join jobs IVAMNs
    join_jobs_ivamns(output_path, max_chunk, noise_cutoff)
