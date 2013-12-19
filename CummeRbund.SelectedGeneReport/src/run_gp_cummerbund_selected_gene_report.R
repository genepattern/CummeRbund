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
suppressPackageStartupMessages(library(tools))
suppressMessages(suppressWarnings(
   library(cummeRbund)
))

option_list <- list(
  make_option("--cuffdiff.input", dest="cuffdiff.input"),
  make_option("--feature.id", dest="feature.id"),
  make_option("--find.similar", dest="find.similar", type="integer", default=NULL),
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
source(file.path(libdir, "gp_cummerbund_selected_gene_report.R"))

report.as.aggregate <- (opts$report.as.aggregate == "yes")
log.transform <- (opts$log.transform == "yes")

check.output.format(opts$output.format)
check.feature.level(opts$feature.level)
genome <- get.genome.from.params(opts$ref.gtf, opts$genome)

if (!is.null(opts$find.similar) && opts$find.similar < 0) {
   stop("If provided, find.similar must be a positive integer")
}

print(c("Running GenePattern CummeRbund Selected Gene Report on data from:", opts$cuffdiff.input))

# Move the following to the util file.
# We 
# Probably need to simply generate the genome arg based on the GTF file name unless the user provides one:
# if (is.null(genome) || grepl('^\\s*$', genome))
#    genome <- sub('\\..+$', '', basename(opts$ref.gtf))
# } 

# Create the job.builder function for run.job
job.builder <- function(cuffdiff.job) {
   GP.CummeRbund.SelectedGene.Report(cuffdiff.job, opts$feature.id, opts$find.similar, opts$ref.gtf, genome,
                                     opts$output.format, opts$feature.level, report.as.aggregate, log.transform)
}

suppressMessages(suppressWarnings(
   run.job(opts$cuffdiff.input, job.builder)
))
