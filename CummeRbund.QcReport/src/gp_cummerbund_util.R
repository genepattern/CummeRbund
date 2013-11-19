## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

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
   stop("Unhandled extension")
}

get.feature.selector <- function(feature.level) {
   # Return a feature selection function based on the specified feature level
   if (feature.level == "gene") { return(genes) }
   if (feature.level == "isoform") { return(isoforms) }
   if (feature.level == "TSS") { return(TSS) }
   if (feature.level == "CDS") { return(CDS) }
}

print.plotObject <- function(plotObj, filename_base, device.open) {
   device.open(filename_base)
   print(plotObj)
   dev.off()  
}