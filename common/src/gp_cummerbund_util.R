## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

check.output.format <- function(output.format) {
   if (!(output.format %in% c("pdf", "svg", "png"))) {
      stop(paste0("Unrecognized output format ", output.format))
   }
}

check.feature.level <- function(feature.level) {
   if (!(feature.level %in% c("genes", "isoforms", "TSS", "CDS"))) {
      stop(paste0("Unrecognized feature level ", feature.level))
   }
}

get.device.open <- function(extension) {
   if (extension == "pdf") {
      return(function(filename_base) {
         pdf(paste0(filename_base, ".pdf"))
      })
   }
   if (extension == "svg") {
      return(function(filename_base) {
         svg(paste0(filename_base, ".svg"))
      })
   }
   if (extension == "png") {
      return(function(filename_base) {
         png(paste0(filename_base, ".png"))
      })
   }
   stop("Unhandled plot file format")
}

get.feature.selector <- function(feature.level) {
   # Return a feature selection function based on the specified feature level
   if (feature.level == "genes") { return(genes) }
   if (feature.level == "isoforms") { return(isoforms) }
   if (feature.level == "TSS") { return(TSS) }
   if (feature.level == "CDS") { return(CDS) }
   stop(paste0("Unrecognized feature level ", feature.level))
}

print.plotObject <- function(plotObj, filename_base, device.open) {
   device.open(filename_base)
   print(plotObj)
   dev.off()  
}

readCufflinks.silent <-function(cuffdiff.job, gtf.file, genome.file) {
   # readCufflinks() is very chatty on stderr even with verbose=FALSE (seems to be due to
   # the underlying RSQLite calls).  Make it quiet...   
   tryCatch({
      sink(file=stdout(), type="message")
      return(readCufflinks(cuffdiff.job, gtfFile = gtf.file, genome = genome.file, verbose=FALSE))
   },
   finally = {
      sink()
   })
}