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

run.job <- function(cuffdiff.input, job.builder) {
   if (file_test("-f", cuffdiff.input)) {
      # If the input is a file, assume it is a SQLite database when passing it to readCufflinks.
      # If it is not, we'll let readCufflinks report the problem.

      # need to create a local symlink.  We don't want to clean this up afterward.
      file.symlink(cuffdiff.input, "cuffData.db")
      job.builder(getwd())
   } 
   else {
      # Otherwise it's a directory.  Check whether it contains a cuffData.db file from a previous
      # CummeRbund job; if so, use that directly.  If not, we'll let readCufflinks process it. 
      dir.contents <- list.files(cuffdiff.input)
      if ("cuffData.db" %in% dir.contents) {
         dbFile <- file.path(cuffdiff.input, "cuffData.db")
         
         # need to create a local symlink.  As above, we don't want to clean this up afterward.
         file.symlink(dbFile, "cuffData.db")
         job.builder(getwd())
      }
      else {
         # We need to set up a local copy of these files (using symlinks) so that the resulting
         # cuffData.db file is created here rather than elsewhere.
         input.job <- "inputJob"
         dir.create(input.job)
         symlinker <- function(x) {
            file.symlink(file.path(cuffdiff.input, x), file.path(input.job, x))
         }
         lapply(dir.contents, FUN = symlinker)

         tryCatch({
            job.builder(input.job)
         },
         finally = {
            # Need to move the SQLite DB to the jobResults dir and clean up the symlinks
            dbFile <- file.path(input.job, "cuffData.db")
            if (file.exists(dbFile)) {
               file.rename(dbFile, file.path(getwd(), "cuffData.db"))
            }
            unlink(input.job, recursive=TRUE)
         })
      }
   }
}

readCufflinks.silent <-function(cuffdiff.job, gtf.file, genome.file) {
   # readCufflinks() is very chatty on stderr even with verbose=FALSE (seems to be due to
   # the underlying RSQLite calls).  Make it quiet...   
   tryCatch({
      sink(file=stdout(), type="message")
      cuff <- readCufflinks(cuffdiff.job, gtfFile = gtf.file, genome = genome.file, verbose=FALSE)
      checkCuffVersionAbove2(cuff)
      print(cuff)
      return(cuff)
   },
   finally = {
      sink()
   })
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
    else {
       print(paste0("  ---> The data was generated with Cufflinks ", myVersionInfo, "."))
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

# Wrapper functions to reduce boilerplate code in the reports
build.XYAxisPlotter <- function(plotTypeName, plotterFunction) {
   function(selected.features, x, y, device.open, filename_base, log.transform=TRUE) {
      plotname <- paste0(filename_base,".",plotTypeName,".",x,"_",y)
      tryCatch({
         plotObj<-plotterFunction(selected.features, x, y, log.transform=TRUE)
         print.plotObject(plotObj, plotname, device.open)
      },
      error = function(err) {
         print(paste0("Error printing the ", plotname, " plot - skipping"))
         print(err)
      })
   }
}

build.standardPlotter <- function(plotTypeName, plotterFunction) {
   function(selected.features, device.open, filename_base, use.replicates=FALSE, log.transform=TRUE) {
      plotname <- paste0(filename_base,".",plotTypeName)
      tryCatch({
         plotObj<-plotterFunction(selected.features, use.replicates=FALSE, log.transform)
         print.plotObject(plotObj, plotname, device.open)
      },
      error = function(err) {
         print(paste0("Error printing the ", plotname, " plot - skipping"))
         print(err)
      })
   }
}

# Utility plot functions for plots used across multiple reports
print.volcanoPlot <- build.XYAxisPlotter("Volcano",
   function(selected.features, x, y, log.transform=TRUE) {
      return(csVolcano(selected.features, x, y, showSignificant=TRUE))
   }
)

print.scatterPlot <- build.XYAxisPlotter("Scatter",
   function(selected.features, x, y, log.transform) {
      return(csScatter(selected.features, x, y, colorByStatus=TRUE, smooth=TRUE, logMode=log.transform))
   }
)


print.MAplot <- build.XYAxisPlotter("MAplot",
   function(selected.features, x, y, log.transform) {
      return(MAplot(selected.features, x, y, smooth=TRUE, logMode=log.transform))
   }
)

print.expressonBarplot <- build.standardPlotter("ExpressionBarplot", 
   function(selected.features, use.replicates, log.transform) {
      return(expressionBarplot(selected.features, replicates=use.replicates, logMode=log.transform))
   }
)

print.expressonPlot <- build.standardPlotter("ExpressionPlot", 
   function(selected.features, use.replicates, log.transform) {
      return(expressionPlot(selected.features, replicates=use.replicates, logMode=log.transform))
   }
)

print.dendrogram <- function(selected.features, device.open, filename_base, use.replicates, log.transform) {
   plotname <- paste0(filename_base,".Dendrogram")
   tryCatch({
      # The dendrogram plots behave differently - apparently the print/plot call is embedded within.
      device.open(plotname)
      dend <- csDendro(selected.features, replicates=use.replicates, logMode=log.transform)
      dev.off()
   },
   error = function(err) {
      print(paste0("Error printing the ", plotname, " plot - skipping"))
      print(err)
   })
}
