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

sessionInfo()

source(file.path(libdir, "gp_cummerbund_util.R"))
source(file.path(libdir, "gp_cummerbund_qc_report.R"))

check.output.format(opts$output.format)
check.feature.level(opts$feature.level)

show.replicates <- (opts$show.replicates == "yes")
log.transform <- (opts$log.transform == "yes")

print(c("Running GenePattern CummeRbund QC Report on data from:", opts$cuffdiff.input))

# Create the job.builder function for run.job
job.builder <- function(cuffdiff.job) {
   GP.CummeRbund.QC.Report(cuffdiff.job, opts$ref.gtf, opts$genome.file, opts$output.format,
                           opts$feature.level, show.replicates, log.transform)
}

suppressMessages(suppressWarnings(
   run.job(opts$cuffdiff.input, job.builder)
))
