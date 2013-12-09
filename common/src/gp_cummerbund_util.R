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
      stop(paste0("Unrecognized output format '", output.format, "'"))
   }
}

check.feature.level <- function(feature.level) {
   if (!(feature.level %in% c("genes", "isoforms", "TSS", "CDS"))) {
      stop(paste0("Unrecognized feature level '", feature.level, "'"))
   }
}

checkCuffVersionAbove2 <- function(cuff) {
  # Emit a warning if the results were generated with a Cufflinks version less than 2.0
  tryCatch({
    print("Checking the version of Cufflinks used to generate this data...")
    myVersionInfo <- runInfo(cuff)[2,2]
    cuffMajorVersion <- as.integer(substring(myVersionInfo, 1, 1))
    if (cuffMajorVersion < 2) {
       print(paste0("The data was generated with Cufflinks ", myVersionInfo, "; Cufflinks 2+ is recommended.  Some plots may not be available."))
    }
  },
  error = function(err) {
    print("Unable to determine the Cufflinks version; Cufflinks 2+ is recommended.  Some plots may not be available.")
  })
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
   stop(paste0("Unhandled plot file format '", extension, "'"))
}

get.feature.selector <- function(feature.level) {
   # Return a feature selection function based on the specified feature level
   if (feature.level == "genes") { return(genes) }
   if (feature.level == "isoforms") { return(isoforms) }
   if (feature.level == "TSS") { return(TSS) }
   if (feature.level == "CDS") { return(CDS) }
   stop(paste0("Unrecognized feature level '", feature.level, "'"))
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

print.volcanoPlot <- function(selected.features, x, y, device.open, filename_base) {
   plotname <- paste0(filename_base,".",x,"_",y)
   tryCatch({
      v<-csVolcano(selected.features, x, y, showSignificant=TRUE)
      print.plotObject(v, plotname, device.open)
   },
   error = function(err) {
      print(paste0("Error printing the ", plotname, " plot - skipping"))
      print(err)
   })
}

print.scatterPlot <- function(selected.features, x, y, log.transform, device.open, filename_base) {
   plotname <- paste0(filename_base,".",x,"_",y)
   tryCatch({
      s <- csScatter(selected.features, x=x, y=y, colorByStatus=TRUE, smooth=TRUE, logMode=log.transform)
      print.plotObject(s, plotname, device.open)
   },
   error = function(err) {
      print(paste0("Error printing the ", plotname, " plot - skipping"))
      print(err)
   })
}

print.MAplot <- function(selected.features, x, y, log.transform, device.open, filename_base) {
   plotname <- paste0(filename_base,".",x,"_",y)
   tryCatch({
      m <- MAplot(selected.features, x=x, y=y, smooth=TRUE, logMode=log.transform)
      print.plotObject(m, plotname, device.open)
   },
   error = function(err) {
      print(paste0("Error printing the ", plotname, " plot - skipping"))
      print(err)
   })
}

print.expressonBarplot <- function(selected.features, show.replicates, log.transform, device.open, filename_base) {
   tryCatch({
      expBarplot<-expressionBarplot(selected.features, replicates=show.replicates, logMode=log.transform)
      print.plotObject(expBarplot, filename_base, device.open)
   },
   error = function(err) {
      print(paste0("Error printing the ", filename_base, " plot - skipping"))
      print(err)
   })
}


print.expressonPlot <- function(selected.features, show.replicates, log.transform, device.open, filename_base) {
   tryCatch({
      expBarplot<-expressionPlot(selected.features, replicates=show.replicates, logMode=log.transform)
      print.plotObject(expPlot, filename_base, device.open)
   },
   error = function(err) {
      print(paste0("Error printing the ", filename_base, " plot - skipping"))
      print(err)
   })
}

print.dendrogram <- function(selected.features, show.replicates, log.transform, device.open, filename_base) {
   tryCatch({
      # The dendrogram plots behave differently - apparently the print/plot call is embedded within.
      device.open(filename_base)
      dend <- csDendro(selected.features, replicates=show.replicates, logMode=log.transform)
      dev.off()
   },
   error = function(err) {
      print(paste0("Error printing the ", filename_base, " plot - skipping"))
      print(err)
   })
}