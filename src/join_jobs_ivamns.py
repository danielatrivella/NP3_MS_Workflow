# script to join the ivamns from different jobs being joined
# first the original msclusterIDs of each IVAMN is mapped to their msclusterID + _ + JOB_CODE
# then the original msclusterIDs with job code will be replaces by the new msclusterID of the joined job
# using the joinedJobsIDs of the clean table
# duplicated edges are merged at the end
# TODO avaliate to run protonated again here
# TODO avaliate how to write the joined annotations to the clean table as columns similar to the annotation script

import os, sys
import pandas as pd
import numpy as np

#  Root Mean Square Error (RMSE) is a metric used to evaluate the difference between predicted and actual values
# used to compute the rt error between annotated nodes
def rmse(predicted, actual):
    return np.sqrt(np.mean((predicted - actual)**2))


# read the join original jobs metadata
# for each original job being joined, read its IVAMN and map its msclusterID to the joined msclusterID
# then concatenate the IVAMNS and merge duplicated edges
def join_jobs_ivamns(output_path):
    if not os.path.isdir(output_path):
        sys.exit("ERROR. The provided output path directory '"+output_path+"' does not exists.")
    output_name = os.path.basename(output_path)

    # read the clean quantifications of the output_name job
    clean_counts_path = os.path.join(output_path, "count_tables", "clean", output_name+"_peak_area_clean.csv")
    if not os.path.isfile(clean_counts_path):
        sys.exit("ERROR. The clean quantification of the current joining jobs '" + clean_counts_path +
                 "' does not exists.")
    clean_counts = pd.read_csv(clean_counts_path, low_memory=False, usecols=['msclusterID', 'joinedJobsIDs',
                                                                             'rtMean', 'rtMin', 'rtMax'])

    # read the default join original jobs metadata, from it extract the Job name, code and path
    original_jobs_metadata_path = os.path.join(output_path, "../..", "join_original_jobs_METADATA.csv")
    if not os.path.isfile(original_jobs_metadata_path):
        sys.exit("ERROR. The metadata table with the original joined jobs '" + original_jobs_metadata_path +
                 "' does not exists.")
    jobs_metadata = pd.read_csv(original_jobs_metadata_path)

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
        else:
            job_code = jobs_metadata.JOB_CODE[i]
            # read the ivamn attribute table
            job_ivamn_att_path = os.path.join(jobs_metadata.JOB_PATH[i], "molecular_networking",
                                                 jobs_metadata.JOBNAME[i]+"_ivamn_attributes.csv")
        if not os.path.isfile(job_ivamn_att_path):
            sys.exit("ERROR. The IVAMN attribute table of the original job with code '"+job_code +
                     "' does not exists. The following path extracted from the joined jobs metadata is not valid: '" +
                     job_ivamn_att_path + "'.")
        job_ivamn_att = pd.read_csv(job_ivamn_att_path, low_memory=False,
                                    usecols=['msclusterID', 'multicharge_ion', 'isotope_ion'])
        if i < jobs_metadata.shape[0]:
            # if this is not the final integration step
            # map the msclusterID_job to the new joinedJobs msclusterID
            job_ivamn_att["msclusterID_job"] = job_ivamn_att.msclusterID.astype(str)+ "_" + job_code
            try:
                job_ivamn_att["msclusterID_new"] = [clean_counts.msclusterID.values[
                     clean_counts.joinedJobsIDs.str.contains("(^"+id + "$|^" + id + ";|;" + id + "$|;" + id + ";)", regex=True)][0]
                 for id in job_ivamn_att.msclusterID_job.values]
            except IndexError:
                # error if any id is not found and indexing the first position of values return an error
                missing_ids = np.asarray([clean_counts.joinedJobsIDs.str.contains("(^" + id + "$|^" + id + ";|;" + id + "$|;" + id + ";)", regex=True).any()
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
        # for each duplicated node, aggregate the multicharge_ion and the isotope_ion columns
        # make them receive the max values among the joined ids
        for dup_node in first_duplicated_nodes:
            dup_id = job_ivamn_att.loc[dup_node, "msclusterID"]
            select_duplicated_rows = (job_ivamn_att.loc[duplicated_joined_nodes, "msclusterID"] == dup_id)
            job_ivamn_att.loc[dup_node, ["multicharge_ion", "isotope_ion"]] = job_ivamn_att.loc[duplicated_joined_nodes,:].loc[
                select_duplicated_rows, ["multicharge_ion", "isotope_ion"]].max()
        #
        # remove the duplicated nodes, except the first duplicated ones
        job_ivamn_att = job_ivamn_att.drop(np.setdiff1d(duplicated_joined_nodes, first_duplicated_nodes), 0)
        # order columns and store the att table
        job_ivamn_att = job_ivamn_att[['msclusterID', 'multicharge_ion', 'isotope_ion']]
        # check if all the joined msclusterIDs are in the clean counts table ids, trought error if not
        if not job_ivamn_att.msclusterID.isin(clean_counts.msclusterID).all():
            sys.exit("ERROR. The joined IVAMN of the original job with code '" + job_code +
                     "' was not correctly created. There is a final msclusterID that do not appear in the joined clean " +
                     " quantification table (problematic IDs: " +
                     ','.join(job_ivamn_att.msclusterID[~job_ivamn_att.msclusterID.isin(clean_counts.msclusterID)]) +
                     "). Something went wrong when joining the IVAMN and reducing redundancies.")
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
                             "). Something went wrong when integrating the joined IVAMNs and reducing redundancies.")
                job_ivamn_att.to_csv(
                    os.path.join(output_path, "molecular_networking", output_name + "_ivamn_attributes_tmp.csv"),
                    index=False)
        # rm current ivamn att table
        del job_ivamn_att

    # TODO at the end recompute rtError and retrieve the new cosine for each available edge


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # print(sys.argv)
        output_path = sys.argv[1]
    else:
        print("Error: One argument must be supplied to join the original IVAMNs of jobs being joined inside the join_jobs command flow:\n",
        " 1 - output_path: Path to the final output directory of the job currently being joined, ",
        "named with the output_name inside the 'outs' directory.\n")
        sys.exit(1)
    # call join jobs IVAMNs
    join_jobs_ivamns(output_path)