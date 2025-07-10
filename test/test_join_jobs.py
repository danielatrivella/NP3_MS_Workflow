# read original clean tables and final joined clean table
# check if all msclusterID from orignal jobs are present in the joinedJobID column of the final joined job

import pandas as pd
import os, sys
import numpy as np

# check the final joined job present in the provides output_path folder, inside the outs dir
# use the original jobs metadata to retrieve the original clean tables, check one job at a time
# and accumulate errors
# check if all original sample codes are present in the final clean table
# check if msclusterIDs and m/zs are maintained (for spectra with basePeakInt above the current min when noise_cutoff was used)
# check if multicharge and isotope ions are maintained
def check_joined_jobs(output_path, noise_cutoff, mz_tolerance=0.025):
    if not os.path.isdir(output_path):
        sys.exit("ERROR. The provided output path directory '" + output_path + "' does not exists.")
    output_name = os.path.basename(output_path)

    print("* Testing the consistency of the Joined Job ", output_name," * ")
    # read the clean quantifications of the output_name job
    clean_counts_path = os.path.join(output_path, "count_tables", "clean", output_name + "_peak_area_clean_ann.csv")
    if not os.path.isfile(clean_counts_path):
        sys.exit("ERROR. The clean quantification of the current joining jobs in terms of peak area with annotations '" +
                 clean_counts_path + "' does not exists.")
    clean_counts = pd.read_csv(clean_counts_path, low_memory=False)

    if noise_cutoff:
        min_basePeakInt = clean_counts.basePeakInt.min()

    # read the default original joined jobs metadata, from it extract the Job name, code and path
    original_jobs_metadata_path = os.path.join(output_path, "../..", "original_jobs_METADATA_JOIN.csv")
    if not os.path.isfile(original_jobs_metadata_path):
        sys.exit("ERROR. The metadata table with the original joined jobs '" + original_jobs_metadata_path +
                 "' does not exists.")
    jobs_metadata = pd.read_csv(original_jobs_metadata_path)
    num_missing_id_job = np.asarray([0]*jobs_metadata.shape[0])
    n_inconsistency = 0

    # for each orignal job, check if all the msclusterIDs are kept in the final joined job
    # check if all m/zs are kept within tolerance (ignore rt) - missing m/z
    for i in range(jobs_metadata.shape[0]):
        job_code = jobs_metadata.JOB_CODE[i]
        print("    - Checking the msclusterIDs, m/z's, multicharge and isotope ions of job code " + job_code)
        # read the clean table
        job_clean_path = os.path.join(jobs_metadata.JOB_PATH[i], "count_tables", "clean",
                                          jobs_metadata.JOBNAME[i] + "_peak_area_clean_ann.csv")

        if not os.path.isfile(job_clean_path):
            sys.exit("ERROR. The clean count table of the original job with code '" + job_code +
                     "' does not exists. The following path extracted from the joined jobs metadata is not valid: '" +
                     job_clean_path + "'.")
        job_clean_counts = pd.read_csv(job_clean_path, low_memory=False,
                                       usecols=['msclusterID', 'mzConsensus', 'rtMean', 'rtMin', 'rtMax',
                                                'multicharge_ion', 'isotope_ion', 'basePeakInt'])
        # remove here the msclusterID with low basePeakInt when noise_cutoff is enabled
        if noise_cutoff:
            job_clean_counts = job_clean_counts.loc[job_clean_counts.basePeakInt >= min_basePeakInt,:]
            # if no spectra is left, skip to next sample
            if job_clean_counts.shape[0] == 0:
                continue

        job_clean_counts["msclusterID_job"] = job_clean_counts.msclusterID.astype(str) + "_" + job_code
        # check if all msclusterIDs were kept
        # search for partial matches with the original IDs
        kept_msclsuterIDs = np.asarray([clean_counts.joinedJobsIDs.str.contains(
                                                    "^" + id + "$|^" + id + ";|;" + id + "$|;" + id + ";", regex=True).any()
                                            for id in job_clean_counts.msclusterID_job.values])
        if not kept_msclsuterIDs.all():
            num_missing_id_job[i] = sum(~kept_msclsuterIDs)
            print("      - ERROR: Missing IDs for job code ",  job_code)
            print("        - A total of ", str(num_missing_id_job[i]),
                  " msclusterIDs were not found in the final joined clean counts, which are: ",
                  ",".join(job_clean_counts.msclusterID_job[~kept_msclsuterIDs].values.astype(str)))

        # check if all kept msclusterIDs maintained their multicharge and isotope ion column values >= than the original
        # 0 -> 1 or 0
        # 1 -> 1
        where_msclsuterIDs = np.concatenate([np.where(clean_counts.joinedJobsIDs.str.contains(
            "^" + id + "$|^" + id + ";|;" + id + "$|;" + id + ";", regex=True))[0]
                                             for id in job_clean_counts.msclusterID_job.values[kept_msclsuterIDs]])
        kept_multicharge_ion = (job_clean_counts.multicharge_ion[kept_msclsuterIDs].values <=
                                clean_counts.multicharge_ion[where_msclsuterIDs].values)
        if not kept_multicharge_ion.all():
            n_inconsistency += sum(~kept_multicharge_ion)
            print("      - ERROR: A total of ", str(sum(~kept_multicharge_ion)),
                  " msclusterIDs with multicharge ions set to 1 from original job code ", job_code,
                  "were not kept as multicharge ions in the final joined counts, which are: ",
                  ",".join(job_clean_counts.msclusterID_job[kept_msclsuterIDs][~kept_multicharge_ion].values.astype(str)))
        kept_isotope_ion = (job_clean_counts.isotope_ion[kept_msclsuterIDs].values <=
                            clean_counts.isotope_ion[where_msclsuterIDs].values)
        if not kept_isotope_ion.all():
            n_inconsistency += sum(~kept_isotope_ion)
            print("      - ERROR: A total of ", str(sum(~kept_isotope_ion)),
                  " msclusterIDs with isotope ions set to 1 from original job code ", job_code,
                  "were not kept as isotope ions in the final joined counts, which are: ",
                  ",".join(job_clean_counts.msclusterID_job[kept_msclsuterIDs][~kept_isotope_ion].values.astype(str)))

        # search for each m/z within tolerance
        kept_mzs = np.asarray([((clean_counts.mzConsensus - 3*mz_tolerance <= mz) & (clean_counts.mzConsensus + mz_tolerance*3 >= mz)).any()
                               for mz in job_clean_counts.mzConsensus])
        if not kept_mzs.all():
            print("      - ERROR: A total of ", str(sum(~kept_mzs)),
                  " m/z consensus from original job code ", job_code,
                  "are missing in the final joined counts, which are: ",
                  ",".join(job_clean_counts.mzConsensus[~kept_mzs].values.astype(str)))
            # sum misiing mzs
            n_inconsistency += sum(~kept_mzs)

    # sum missing mscluster ids
    n_inconsistency += sum(num_missing_id_job)
    print("    - Checking the original sample codes ")
    # check if all sample code is present in the final joined clean table
    # read the default original joined jobs samples metadata, from it extract the original samples codes
    original_samples_metadata_path = os.path.join(output_path, "../..", "original_samples_METADATA.csv")
    if not os.path.isfile(original_samples_metadata_path):
        sys.exit("ERROR. The metadata table with the original samples of the joined jobs '" + original_samples_metadata_path +
                 "' does not exists.")
    samples_metadata = pd.read_csv(original_samples_metadata_path)

    count_column_names = clean_counts.columns[clean_counts.columns.str.endswith("_area")].values
    samples_metadata["SAMPLE_CODE_AREA"] = samples_metadata.SAMPLE_CODE.astype(str) + "_area"
    kept_sample_code = samples_metadata.SAMPLE_CODE_AREA.isin(count_column_names)
    if not kept_sample_code.all():
        print("      - ERROR: A total of ", str(sum(~kept_sample_code)),
              " original sample codes are missing, which are: ",
              ",".join(samples_metadata.SAMPLE_CODE[~kept_sample_code]))

    n_inconsistency += sum(~kept_sample_code)
    if n_inconsistency == 0:
        print("Done! :)\n")
    else:
        print("ERROR! A total of", str(n_inconsistency), "inconsistencies were found. :(\n")


if __name__ == "__main__":
    mz_tol = 0.025
    if len(sys.argv) > 2:
        # print(sys.argv)
        output_path = sys.argv[1]
        noise_cutoff = bool(float(sys.argv[2]))
        if len(sys.argv) > 3:
            mz_tol = float(sys.argv[3])
    else:
        print("Error: Two arguments must be supplied to check the final joined jobs consistency within the original jobs ids and sample codes:\n",
              " 1 - output_path: Path to the final output directory of the joined job, ",
              "named with the output_name inside the 'outs' directory;\n",
              " 2 - noise_cutoff: a boolean True or False indicating if the noise_cutoff was used in the np3 join_jobs processing "
              "and thus spectra with a low basePeakInt (values smaller the current minimum basePeakInt) should be "
              "ignored from the original jobs being checked (spectra that were probable removed by the noise cutoff).\n"
              " 3 - mz_tolerance: the tolerance in Daltons to assume two precursor masses the same (default: 0.025).")
        sys.exit(1)
    # call check joined job
    check_joined_jobs(output_path, noise_cutoff, mz_tol)
