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

suppressMessages(suppressWarnings(
   library(optparse))
))
suppressMessages(suppressWarnings(
   library(cummeRbund)
))

# Based on info from Loyal Goff, the ref.gtf and genome parameters are unused at this time.
# There is reason to believe that these might be brought back, however, so they have only
# been removed from the manifest and not from the R code.
option_list <- list(
  make_option("--cuffdiff.input", dest="cuffdiff.input"),
  make_option("--geneset.file", dest="geneset.file"),
  make_option("--selected.conditions", dest="selected.conditions", default=NULL),
  make_option("--ref.gtf", dest="ref.gtf", default=NULL),
  make_option("--genome", dest="genome", default=NULL),
  make_option("--output.format", dest="output.format"),
  make_option("--feature.level", dest="feature.level"),
  make_option("--report.as.aggregate", dest="report.as.aggregate"),
  make_option("--log.transform", dest="log.transform")
  )

opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=args)
print(opt)
opts <- opt$options

sessionInfo()

source(file.path(libdir, "gp_cummerbund_util.R"))
source(file.path(libdir, "gp_cummerbund_geneset_report.R"))

report.as.aggregate <- (opts$report.as.aggregate == "yes")
log.transform <- (opts$log.transform == "yes")

check.output.format(opts$output.format)
check.feature.level(opts$feature.level)
genome <- get.genome.from.params(opts$ref.gtf, opts$genome)

selected.conditions <- get.selected.conditions(opts$selected.conditions)

print(c("Running GenePattern CummeRbund Geneset Report with data from:", opts$cuffdiff.input))

# Create the job.builder function for run.job
job.builder <- function(cuffdiff.job) {
   GP.CummeRbund.Geneset.Report(cuffdiff.job, opts$geneset.file, selected.conditions, opts$ref.gtf, genome,
                                opts$output.format, opts$feature.level, report.as.aggregate, log.transform)
}

suppressMessages(suppressWarnings(
   run.job(opts$cuffdiff.input, job.builder)
))