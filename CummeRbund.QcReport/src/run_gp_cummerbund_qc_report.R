## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

args <- commandArgs(trailingOnly=TRUE)

libdir <- args[1]
site.library <- args[2]

cat("\nLibrary dir: ",site.library)
.libPaths(site.library)

suppressPackageStartupMessages(library(optparse))
suppressMessages(suppressWarnings(
   library(cummeRbund)
))

option_list <- list(
  make_option("--cuffdiff.input", dest="cuffdiff.input"),
  make_option("--ref.gtf", dest="ref.gtf", default=NULL),
  make_option("--genome.file", dest="genome.file", default=NULL),
  make_option("--output.format", dest="output.format"),
  make_option("--feature.level", dest="feature.level"),
  make_option("--show.replicates", dest="show.replicates"),
  make_option("--log.transform", dest="log.transform")
  )

opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=args)
print(opt)
opts <- opt$options

if (!(opts$output.format %in% c("pdf", "svg", "png"))) {
   stop(paste0("Unrecognized output format ", opts$output.format))
}

if (!(opts$feature.level %in% c("gene", "isoform", "TSS", "CDS"))) {
   stop(paste0("Unrecognized feature level ", opts$feature.level))
}

show.replicates <- (opts$show.replicates == "yes")
log.transform <- (opts$log.transform == "yes")

source(file.path(libdir, "gp_cummerbund_util.R"))
source(file.path(libdir, "gp_cummerbund_qc_report.R"))

sessionInfo()

run.job <- function(cuffdiff.input, gtf.file, genome.file, output.format, feature.level, show.replicates, log.transform) {
   print(c("Running GenePattern CummeRbund QC Report on data from:", basename(cuffdiff.input)))

   if (file_test("-f", cuffdiff.input)) {
      # If the input is a file, assume it is a SQLite database when passing it to readCufflinks.
      # If it is not, we'll let readCufflinks report the problem.

      # need to create a local symlink.  We don't want to clean this up afterward.
      file.symlink(cuffdiff.input, "cuffData.db")
      
      GP.CummeRbund.QC.Report(cuffdiff.job = getwd(), gtf.file, genome.file, output.format, 
                              feature.level, show.replicates, log.transform)
   } 
   else {
      # Otherwise it's a directory.  Check whether it contains a cuffData.db file from a previous
      # CummeRbund job; if so, use that directly.  If not, we'll let readCufflinks process it. 
      dir.contents <- list.files(cuffdiff.input)
      if ("cuffData.db" %in% dir.contents) {
         dbFile <- file.path(cuffdiff.input, "cuffData.db")
         
         # need to create a local symlink.  As above, we don't want to clean this up afterward.
         file.symlink(dbFile, "cuffData.db")
         
         GP.CummeRbund.QC.Report(cuffdiff.job = getwd(), gtf.file, genome.file, output.format,
                                 feature.level, show.replicates, log.transform)
      }
      else {
         # We need to set up a local copy of these files (using symlinks) so that the resulting
         # cuffData.db file is created here rather than elsewhere.
         input.job <- "inputJob"
         dir.create(input.job)
         symlinker <- function(x) {
            file.symlink(file.path(cuffdiff.input, x), file.path(input.job, x))
         }
         apply(dir.contents, FUN = symlinker)

         tryCatch({
            GP.CummeRbund.QC.Report(cuffdiff.job = input.job, gtf.file, genome.file, output.format,
                                    feature.level, show.replicates, log.transform)
         },
         finally = {
            # Need to move the SQLite DB to the jobResults dir and clean up the symlinks
            dbFile <- file.path(input.job, "cuffData.db")
            if (file.exists(dbFile)) {
               file.rename(dbFile, file.path(getwd(), "cuffData.db"))
            }
            unlink(input.job)
         })
      }
   }
}

suppressMessages(suppressWarnings(
   run.job(opts$cuffdiff.input, opts$gtf.file, opts$genome.file, opts$output.format,
           opts$feature.level, show.replicates, log.transform)
))
