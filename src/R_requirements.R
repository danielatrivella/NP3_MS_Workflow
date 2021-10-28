depend_pcks <- c('purrr', 'dplyr', 'stats', 'inline', 'readr',  # 'rJava',
                 'RcppArmadillo', "gtools", 'BiocManager', 'stringr')
depend_pcks_version <- c("0.3.2", "0.8.3", "3.6.3", "0.3.16", "1.3.1", # "0.9-13",
                         "0.9.700.2.0", "3.8.2", "1.30.10", "1.4.0")
bio_pcks <- c('Rcpp','mzR', 'MSnbase', 'xcms')

# install devtools first to then install the packages in the desired versions
if (!require(devtools)) {
  cat('\n* Installing R package: devtools *\n')
  install.packages("devtools", repos = 'https://cloud.r-project.org/', 
                   dependencies = TRUE)
}

# install the depend packages using the devtools install version
depend_pcks_missing <- setdiff( depend_pcks, rownames(installed.packages()))
if (length(depend_pcks_missing) > 0)  
{ 
  depend_pcks_version <- depend_pcks_version[match(depend_pcks_missing, depend_pcks)]
  cat('\n* Installing R packages: ', paste(depend_pcks_missing, collapse = ", "), 
      '* \n')
  for (i in 1:length(depend_pcks_missing)) {
    install_version(depend_pcks_missing[i], version = depend_pcks_version[i],
                    repos = 'https://cloud.r-project.org/', 
                    dependencies = TRUE)
  }
}

# install metfrag
# if (length(setdiff(c("rJava", "metfRag"), rownames(installed.packages()))) > 0 && 
#     all(setdiff(c("rJava", "metfRag"), rownames(installed.packages())) == "metfRag"))
# {
#   require("rJava")
#   require("devtools") 
#   cat('\n* Installing R package: metfRag *\n')
#   tryCatch(install_github("c-ruttkies/MetFragR/metfRag"),
#            error=function(e) print(e))
#   if (!require("metfRag")) 
#   {
#       tryCatch(install.packages("metfRag",repos=NULL,type="source"),
#            error=function(e) print(e))
#   }
# }

# install the bioconductor packages 
bio_pcks_missing <- setdiff(bio_pcks, rownames(installed.packages()))
if (length(bio_pcks_missing) > 0)
{
  # upgrade BiocManager and install dependencies
  if(.Platform$OS.type == "unix") {
	  cat('\n* Installing R packages: \n', paste(bio_pcks_missing, collapse = ", "), 
		  '*\n')
	  BiocManager::install(version="3.10") 
	  tryCatch(BiocManager::install(bio_pcks_missing, 
									update = FALSE, version="3.10"),
			   error=function(e) print(e))
  } else {
	# if on windows, install all bio packeges using biocmanager
	cat('\n* Installing R packages: \n', paste(bio_pcks, collapse = ", "), 
		  '*\n')
	BiocManager::install(version="3.10") 
	tryCatch(BiocManager::install(bio_pcks, 
									update = FALSE, version="3.10"),
			   error=function(e) print(e))
  }
}
# depend_pcks <- c('devtools', depend_pcks, 'metfRag', bio_pcks)
depend_pcks <- c('devtools', depend_pcks, bio_pcks)

# check if the packages were installed
if (length(setdiff(depend_pcks, rownames(installed.packages()))) > 0) 
{
  stop("\n\nThe following R dependencies could not be installed:\n ", 
       paste(setdiff(depend_pcks, rownames(installed.packages())), 
             collapse = ", "),
       "\n\nLook for ERROR messages, install the missing libraries and retry.")
}

# check if the packages can be loaded
cat("\n\n* Checking if the R packages can be loaded *\n\n")
required_pcks <- sapply(depend_pcks, 
                        function(x)
                            tryCatch(require(x ,character.only = TRUE), 
                                     error=function(e) FALSE))
if (any(!required_pcks)) {
  stop("The following R dependencies could not be loaded: \n", 
       paste(depend_pcks[!required_pcks], collapse = ", "),
       "\nLook for ERROR messages, install the missing libraries and retry.")
} else {
  cat('\n  * All R packages could be loaded successfully! *\n')
}