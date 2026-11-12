# final set of selected descriptors - top 24
set_reference_descriptors_bestCos_top24 <- c('WPOL','ATSp1','ATSp2','SP.3','Zagreb','VP.3','naAromAtom','SP.2','ATSp3','VP.2','ATSm3','nB','ATSm2','nAromRings','nAromBlocks','C2SP2','khs.aasC','nAtomP','khs.aaCH','LipinskiFailures','MW','nAtom','VP.1','WTPT.1')
# from MPH; set_descriptors_enamine = c("naAromAtom","WPOL","Zagreb", "nAromBond", "ATSp1", "ATSp2", "nB", "nAromRings", "ATSm3","ATSm4","ATSp3","ATSp5", "nAromBlocks", "ATSm5", "nAtom","VAdjMat","C2SP2","khs.sssCH","MDEC.34","LipinskiFailures")
# from the bacs article:  ALogp2, C1SP2, C1SP3, C3SP2, C3SP3, LipinskiFailures, MDEC.22, MDEC.33, MDEN.12, MDEN.22, MDEN.23, MDEO.22, nAcid, nBase, nRingBlocks, nRings7, WTPT.5 and XLogP
# from the rcdk - all: 
set_descriptors_rcdk <- c('nSmallRings','nAromRings','nRingBlocks','nAromBlocks','nRings3','nRings4','nRings5','nRings6','nRings7','nRings8','nRings9','tpsaEfficiency','Zagreb','WPATH','WPOL','WTPT.1','WTPT.2','WTPT.3','WTPT.4','WTPT.5','VAdjMat','VABC','TopoPSA','topoShape','geomShape','PetitjeanNumber','MDEC.11','MDEC.12','MDEC.13','MDEC.14','MDEC.22','MDEC.23','MDEC.24','MDEC.33','MDEC.34','MDEC.44','MDEO.11','MDEO.12','MDEO.22','MDEN.11','MDEN.12','MDEN.13','MDEN.22','MDEN.23','MDEN.33','khs.sLi','khs.ssBe','khs.ssssBe','khs.ssBH','khs.sssB','khs.ssssB','khs.sCH3','khs.dCH2','khs.ssCH2','khs.tCH','khs.dsCH','khs.aaCH','khs.sssCH','khs.ddC','khs.tsC','khs.dssC','khs.aasC','khs.aaaC','khs.ssssC','khs.sNH3','khs.sNH2','khs.ssNH2','khs.dNH','khs.ssNH','khs.aaNH','khs.tN','khs.sssNH','khs.dsN','khs.aaN','khs.sssN','khs.ddsN','khs.aasN','khs.ssssN','khs.sOH','khs.dO','khs.ssO','khs.aaO','khs.sF','khs.sSiH3','khs.ssSiH2','khs.sssSiH','khs.ssssSi','khs.sPH2','khs.ssPH','khs.sssP','khs.dsssP','khs.sssssP','khs.sSH','khs.dS','khs.ssS','khs.aaS','khs.dssS','khs.ddssS','khs.sCl','khs.sGeH3','khs.ssGeH2','khs.sssGeH','khs.ssssGe','khs.sAsH2','khs.ssAsH','khs.sssAs','khs.sssdAs','khs.sssssAs','khs.sSeH','khs.dSe','khs.ssSe','khs.aaSe','khs.dssSe','khs.ddssSe','khs.sBr','khs.sSnH3','khs.ssSnH2','khs.sssSnH','khs.ssssSn','khs.sI','khs.sPbH3','khs.ssPbH2','khs.sssPbH','khs.ssssPb','Kier1','Kier2','Kier3','HybRatio','fragC','FMF','ECCEN','SP.0','SP.1','SP.2','SP.3','SP.4','SP.5','SP.6','SP.7','VP.0','VP.1','VP.2','VP.3','VP.4','VP.5','VP.6','VP.7','SPC.4','SPC.5','SPC.6','VPC.4','VPC.5','VPC.6','SC.3','SC.4','SC.5','SC.6','VC.3','VC.4','VC.5','VC.6','SCH.3','SCH.4','SCH.5','SCH.6','SCH.7','VCH.3','VCH.4','VCH.5','VCH.6','VCH.7','C1SP1','C2SP1','C1SP2','C2SP2','C3SP2','C1SP3','C2SP3','C3SP3','C4SP3','ATSp1','ATSp2','ATSp3','ATSp4','ATSp5','ATSm1','ATSm2','ATSm3','ATSm4','ATSm5','ATSc1','ATSc2','ATSc3','ATSc4','ATSc5','topoShape1','geomShape1','MOMI.X','MOMI.Y','MOMI.Z','MOMI.XY','MOMI.XZ','MOMI.YZ','MOMI.R','LOBMAX','LOBMIN','GRAV.1','GRAV.2','GRAV.3','GRAVH.1','GRAVH.2','GRAVH.3','GRAV.4','GRAV.5','GRAV.6','PPSA.1','PPSA.2','PPSA.3','PNSA.1','PNSA.2','PNSA.3','DPSA.1','DPSA.2','DPSA.3','FPSA.1','FPSA.2','FPSA.3','FNSA.1','FNSA.2','FNSA.3','WPSA.1','WPSA.2','WPSA.3','WNSA.1','WNSA.2','WNSA.3','RPCG','RNCG','RPCS','RNCS','THSA','TPSA','RHSA','RPSA','Fsp3','XLogP','MW','LipinskiFailures','nRotB','MLogP','nAtomP','nAtomLC','nB','nBase','nAtom','nAromBond','naAromAtom','ALogP','ALogp2','AMR','nAcid')
set_reference_topological_descriptors_names_top24 <- c("org.openscience.cdk.qsar.descriptors.molecular.SmallRingDescriptor",
                                                       "org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor",
                                                       "org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor",
                                                       "org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor",
                                                       "org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor",
                                                       "org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor",
                                                       "org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor",
                                                       "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability",
                                                       "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass")
set_reference_constitutional_descriptors_names_top24 <- c("org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
                                                          "org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor",
                                                          "org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor",
                                                          "org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor",
                                                          "org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor",
                                                          "org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor")
set_reference_geometrical_descriptors_names_top24 <- c()
set_reference_descriptors_names_top24 <- c(set_reference_topological_descriptors_names_top24,
                                     set_reference_constitutional_descriptors_names_top24,
                                     set_reference_geometrical_descriptors_names_top24)

# change global java options memory - set maximum heap size to 2Gb - allow processing big SMILES
options(java.parameters = "-Xmx2g")
# load libraries
library(rJava)
suppressPackageStartupMessages(library(rcdk))
suppressPackageStartupMessages(library(dplyr))
library(readr)

# read input
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Two arguments must be supplied to calculate the RCDK descriptors (topological, geometrical and constitutional) from a list of SMILES (valid only - not NA), the table containing the smiles and the column name. Optionally, the output path may be informed:\n",
       " 1 - smiles_file: the path to the table in CSV format containing SMILES string in a column;\n",
       " 2 - smiles_col: the name of the column containing the SMILES string;\n",
       " 3 - output_path: (optional) the output directory to store the table with the descriptors columns concatenated to the smiles_file",
       call.=FALSE)
} else {
  smiles_file <- file.path(args[[1]])
  smiles_col <- args[[2]]
  parallel_cores <- as.numeric(args[[3]])
  if (parallel_cores < 1) {
    cat("Invalid parallel_cores value (", parallel_cores, 
        "). The number of cores to be used for parallel processing must be at least 2 ",
        "or 1 for disabling parallelization. The number of cores was set to 1 by default.\n")
    warning("Invalid parallel_cores parameter value.")
    parallel_cores <- 1    
  } else if (parallel_cores > 1 && !require(parallel)) {
    cat("Parallel library is not available, disabling parallelization. ",
        "The number of cores was set to 1 by default.\n")
    warning("Invalid parallel_cores parameter value.")
    parallel_cores <- 1 
  } else if (parallel_cores > 1 && !is.na(detectCores(logical = FALSE)) && parallel_cores >= detectCores(logical = FALSE)) {
    cat("The number of paralell_cores =",parallel_cores,
        " is greater or equal than the number of physical cores =",detectCores(logical = FALSE),
        " in your machine. This is not recommended,",
        "the number of parallel clusters will be limited by the number of physical cores minus two =",
        detectCores(logical = FALSE)-2, ".\n")
    warning("Invalid parallel_cores parameter value.")
    parallel_cores <- detectCores(logical = FALSE)-2
  }
  output_path <- ""
  if (length(args) >= 4) {
    output_path <- file.path(args[[4]])
    if (!dir.exists(output_path))
    {
    stop("The path to the output directory '", output_path,
         "' do not exists. Provide a valid path to store the results.")
    }
  }
}


compute_descriptors_smiles <- function(smiles_char, 
                                       topological_desc_names, 
                                       geometrical_desc_names,
                                       constitutional_desc_names,
                                       set_reference_descriptors_cols_list) {
  parsed_smiles <- parse.smiles(smiles_char)
  if (parsed_smiles == "NULL") {
    descriptors_data <- as.list(rep.int(NA_real_, length(set_reference_descriptors_cols_list)),all.names=set_reference_descriptors_cols_list)
    names(descriptors_data) <- set_reference_descriptors_cols_list
    return(descriptors_data)
  }
  #cat("    - Calculating topological descriptors \n")
  # evaluate the descriptors
  topological_descriptors <- tryCatch(eval.desc(parsed_smiles, 
                                                topological_desc_names), 
                                      error = function(e) {
                                        #cat(paste(e, "- SMILES", i, ":", smiles_char_list[i], "\n"))
                                        NULL
                                      })
  if (!is.null(topological_descriptors) && 
      (nrow(topological_descriptors) == 0 || ncol(topological_descriptors) == 0)) {
    topological_descriptors <- NULL
  }
  #cat("    - Calculating geometrical descriptors\n")
  # evaluate the descriptors
  geometrical_descriptors <- tryCatch(eval.desc(parsed_smiles, 
                                                geometrical_desc_names), 
                                      error = function(e) {
                                        #cat(paste(e, "- SMILES", i, ":", smiles_char_list[i], "\n"))
                                        NULL
                                      })
  if (!is.null(geometrical_descriptors) && 
      (nrow(geometrical_descriptors) == 0 || ncol(geometrical_descriptors) == 0)) {
    geometrical_descriptors <- NULL
  }
  #cat("    - Calculating constitutional descriptors\n")
  # evaluate the descriptors
  constitutional_descriptors <- tryCatch(eval.desc(parsed_smiles, 
                                                   constitutional_desc_names), 
                                         error = function(e) {
                                           #cat(paste(e, "- SMILES", i, ":", smiles_char_list[i], "\n"))
                                           NULL
                                         })
  if (!is.null(constitutional_descriptors) && 
      (nrow(constitutional_descriptors) == 0 || ncol(constitutional_descriptors) == 0)) {
    constitutional_descriptors <- NULL
  }
  # concatenate the descriptors result for current smiles
  # if any valid result, otherwise leave its data as NA
  descriptors_data <- bind_cols(topological_descriptors, geometrical_descriptors, constitutional_descriptors)
  if (!is.null(descriptors_data) && nrow(descriptors_data) != 0 && 
      ncol(descriptors_data) != 0) {
    # filter the selected valid descriptors
    descriptors_data <- descriptors_data[,names(descriptors_data) %in% 
                                           set_reference_descriptors_cols_list]
    # and set the selected descriptors to the current smiles in the data table
    #data_SMILES[i,names(descriptors_data)] <- descriptors_data
  } else {
    descriptors_data <- as.list(rep.int(NA_real_, length(set_reference_descriptors_cols_list)),all.names=set_reference_descriptors_cols_list)
    names(descriptors_data) <- set_reference_descriptors_cols_list
  }
  return(descriptors_data)
}

## Extract the rcdk descriptors (topological, geometrical and constitutional) - all descriptors are retrieved by default -
# of valid SMILES (not NA) present in a table, where *smiles_file* is the path 
# to this table in CSV format containing SMILES string in a 
# column named with the *smiles_col* name
# if set_reference_descriptors_names is not 'all' it should contain a list with the
# descriptors names to be calculated;
# if set_reference_descriptors_list is not 'all' it should contain a list with the descriptors
# resulting column names to be kept, otherwise the set_descriptors_rcdk is used which contains all descriptors column names
calculate_rcdk_descriptors <- function(smiles_file, smiles_col="SMILES", set_reference_descriptors_names="all",
                                       set_reference_descriptors_cols_list="all", output_path="", parallel_cores) {
  # check if smiles_file exists and read it
  smiles_file <- file.path(smiles_file)
  if (!file.exists(smiles_file))
  {
    stop("The table '", smiles_file, "' with the list of smiles do not exists. ",
         "Provide a valid path to where the csv file with the list ",
         "of smiles to calculate the RCDK descriptors is located.")
  }
  data <- suppressMessages(read_csv(file=smiles_file,guess_max = 5000))
  # check if the smiles column is present
  if (!(smiles_col %in% names(data))) {
    stop("The provided table '", smiles_file, "' does not contain the smiles ",
         "column named '", smiles_col, 
         "'. Provide a valid name for the column containing the SMILES list.")
  }
  ti <- Sys.time()
  # select the valid smiles
  valid_smiles <- !is.na(data[[smiles_col]])
  n <- sum(valid_smiles)
  cat(paste("\n  - Extracted valid SMILES:",n," entries\n"))
  # parse the valid smiles using rcdk to generate molecules info
  smiles_char <- as.character(data[valid_smiles, smiles_col][[1]])
  
  # initialize the selected descriptors cols with NA values type float
  # if all descriptors were selected, set the reference list to all rcdk descriptors
  if (length(set_reference_descriptors_cols_list) == 1 && set_reference_descriptors_cols_list == "all")
  {
    set_reference_descriptors_cols_list <- set_descriptors_rcdk
  } 
  
  cat(paste("  - Calculating",length(set_reference_descriptors_cols_list),"descriptors \n"))
  cat(paste("    - Number of parallel cores to use:",parallel_cores,"\n"))
  # calculate the descriptors for the valid smiles, parse them in each iteration
  # make this one smiles per round, with try catch errors
  # Topological descriptors
  topological_desc_names <- get.desc.names("topological")
  if (length(set_reference_descriptors_names) > 1 || set_reference_descriptors_names != "all")
  {
    # filter the descriptors names
    topological_desc_names <- topological_desc_names[topological_desc_names %in% set_reference_descriptors_names]
  }
  # Geometrical descriptors
  geometrical_desc_names <- get.desc.names("geometrical")
  if (length(set_reference_descriptors_names) > 1 || set_reference_descriptors_names != "all")
  {
    # filter the descriptors names
    geometrical_desc_names <- geometrical_desc_names[geometrical_desc_names %in% set_reference_descriptors_names]
  }
  # Constitutional descriptors
  constitutional_desc_names <- get.desc.names("constitutional")
  problematic_desc <- "LongestAliphaticChainDescriptor"
  constitutional_desc_names <- constitutional_desc_names[!endsWith(constitutional_desc_names, problematic_desc)]
  if (length(set_reference_descriptors_names) > 1 || set_reference_descriptors_names != "all")
  {
    # filter the descriptors names
    constitutional_desc_names <- constitutional_desc_names[constitutional_desc_names %in% set_reference_descriptors_names]
  }
  if (parallel_cores > 1 & length(smiles_char) > parallel_cores) {
    # using makeCluster here because rJava can not be parallelized with mclapply
    # Note: Do NOT export 'parsed_smiles' because Java objects break across PSOCK workers
    # 1. Create the PSOCK cluster
    cl <- makeCluster(parallel_cores, type = "PSOCK")

    # 2. Load rcdk (and any other packages needed inside compute_descriptors_smiles) on each worker
    clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(rcdk))
      library(dplyr)
    })

    # 3. Export your function and static variables to all workers
    clusterExport(
      cl,
      varlist = c(
        "compute_descriptors_smiles",
        "smiles_char",
        "topological_desc_names",
        "geometrical_desc_names",
        "constitutional_desc_names",
        "set_reference_descriptors_cols_list"
      ),
      envir = environment()
    )

    # 4. Run parLapply and bind the results
    res_data <- bind_rows(
      parLapply(
        cl,
        smiles_char,
        compute_descriptors_smiles,
        topological_desc_names = topological_desc_names,
        geometrical_desc_names = geometrical_desc_names,
        constitutional_desc_names = constitutional_desc_names,
        set_reference_descriptors_cols_list = set_reference_descriptors_cols_list
      )
    )

    # 5. Stop the cluster to free RAM
    stopCluster(cl)
  } else {
    res_data <- bind_rows(lapply(smiles_char, 
                                 compute_descriptors_smiles, 
                                 topological_desc_names=topological_desc_names, 
                                 geometrical_desc_names=geometrical_desc_names,
                                 constitutional_desc_names=constitutional_desc_names,
                                 set_reference_descriptors_cols_list=set_reference_descriptors_cols_list))
  }
  # check size matches
  if (nrow(res_data) != length(smiles_char)) {
    stop("Error when computing the descriptors, wrong number of final computed descriptors list.")
  }
  
  tf <- Sys.time()
  cat(paste("Done computing", ncol(res_data),"RCDK descriptors for the valid SMILES in", 
      round(tf-ti, 2), units(tf-ti), "!\n\n"))
  
  desc_error <- rowSums(is.na(res_data))
  if (any(desc_error > 0)) {
    cat(paste("    - A total of", sum(desc_error > 0), 
              "(",round(sum(desc_error > 0)/nrow(res_data)*100,1),"%) valid SMILES did not have all the descriptors computed successfully.",
              "They are listed below:\n"))
    desc_error <- sapply(which(desc_error > 0), function(i) 
      cat(paste("    - SMILES", i, "(missing",desc_error[i],"descriptors):", smiles_char[i], " \n")))
    
  }
  
  # add the descriptors result to the original data table
  data[valid_smiles, names(res_data)] <- res_data
  # save the original data table with the descriptors columns
  if (output_path == "") {
    write_csv(data, sub(".csv", "_descriptorsCDK.csv", smiles_file))
  } else {
    write_csv(data, file.path(output_path, sub(".csv", "_descriptorsCDK.csv", basename(smiles_file))))
  }
}

# call rcdk computation
calculate_rcdk_descriptors(smiles_file = smiles_file, smiles_col = smiles_col,
                           set_reference_descriptors_names=set_reference_descriptors_names_top24,
                           set_reference_descriptors_cols_list=set_reference_descriptors_bestCos_top24,
                           output_path=output_path, parallel_cores=parallel_cores)