## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2015) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

suppressMessages(suppressWarnings(library(survival)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(colorspace)))
suppressMessages(suppressWarnings(library(dichromat)))
suppressMessages(suppressWarnings(library(munsell)))
suppressMessages(suppressWarnings(library(labeling)))
suppressMessages(suppressWarnings(library(bitops)))
suppressMessages(suppressWarnings(library(Formula)))
suppressMessages(suppressWarnings(library(DBI)))
suppressMessages(suppressWarnings(library(digest)))
suppressMessages(suppressWarnings(library(gtable)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(Rcpp)))
suppressMessages(suppressWarnings(library(plyr)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(lattice)))
suppressMessages(suppressWarnings(library(latticeExtra)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(proto)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(Hmisc)))
suppressMessages(suppressWarnings(library(XML)))
suppressMessages(suppressWarnings(library(RCurl)))
suppressMessages(suppressWarnings(library(RSQLite)))
suppressMessages(suppressWarnings(library(fastcluster)))
suppressMessages(suppressWarnings(library(BiocGenerics)))
suppressMessages(suppressWarnings(library(IRanges)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(BSgenome)))
suppressMessages(suppressWarnings(library(zlibbioc)))
suppressMessages(suppressWarnings(library(Rsamtools)))
suppressMessages(suppressWarnings(library(biomaRt)))
suppressMessages(suppressWarnings(library(Biobase)))
suppressMessages(suppressWarnings(library(AnnotationDbi)))
suppressMessages(suppressWarnings(library(rtracklayer)))
suppressMessages(suppressWarnings(library(GenomicFeatures)))
suppressMessages(suppressWarnings(library(biovizBase)))
suppressMessages(suppressWarnings(library(Gviz)))
suppressMessages(suppressWarnings(library(cummeRbund)))

sessionInfo()

args <- commandArgs(trailingOnly=TRUE)

libdir <- args[1]

# Based on info from Loyal Goff, the ref.gtf and genome parameters are unused at this time.
# There is reason to believe that these might be brought back, however, so they have only
# been removed from the manifest and not from the R code.
option_list <- list(
  make_option("--cuffdiff.input", dest="cuffdiff.input"),
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

source(file.path(libdir, "gp_cummerbund_util.R"))
source(file.path(libdir, "gp_cummerbund_selected_conditions_report.R"))

report.as.aggregate <- (opts$report.as.aggregate == "yes")
log.transform <- (opts$log.transform == "yes")

check.output.format(opts$output.format)
check.feature.level(opts$feature.level)
genome <- get.genome.from.params(opts$ref.gtf, opts$genome)

selected.conditions <- get.selected.conditions(opts$selected.conditions)
if (NROW(selected.conditions) == 0) {
   stop("No conditions specified.  Please specify at least one condition.")
}

print(c("Running GenePattern CummeRbund Selected Conditions Report with data from:", opts$cuffdiff.input))

# Create the job.builder function for run.job
job.builder <- function(cuffdiff.job) {
   GP.CummeRbund.SelectedConditions.Report(cuffdiff.job, selected.conditions, opts$ref.gtf, genome, opts$output.format, 
                                           opts$feature.level, report.as.aggregate, log.transform)
}

suppressMessages(suppressWarnings(
   run.job(opts$cuffdiff.input, job.builder)
))

sessionInfo()