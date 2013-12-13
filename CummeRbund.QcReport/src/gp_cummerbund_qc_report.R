## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GP.CummeRbund.QC.Report <- function(cuffdiff.job, gtf.file, genome.file, output.format,
                                    feature.level, show.replicates, log.transform) {
   device.open <- get.device.open(output.format)
   
   feature.selector <- get.feature.selector(feature.level)

   cuff <- readCufflinks.silent(cuffdiff.job, gtf.file, genome.file)
   selected.features <- feature.selector(cuff)
   
   # Write out a table of the differentially expressed features
   write.significance.data(selected.features, diffData, 
                           paste0('QC.sig_diffExp_', feature.level,'.txt'))

   # Write out tables of the significant promoters, splicing, and relCDS data
   write.significance.data(promoters(cuff), distValues, 'QC.sig_promoter_data.txt')
   write.significance.data(splicing(cuff), distValues, 'QC.sig_splicing_data.txt')
   write.significance.data(relCDS(cuff), distValues, 'QC.sig_relCDS_data.txt')

   # Write out the various plots
   print.dispersionPlot(selected.features, device.open, "QC")
   print.fpkmSCVPlot(selected.features, device.open, "QC")
   print.densityPlot(selected.features, device.open, "QC", show.replicates, log.transform)
   print.Boxplot(selected.features, device.open, "QC", show.replicates, log.transform)
   print.scatterMatrix(selected.features, device.open, "QC", show.replicates, log.transform)
   print.volcanoMatrix(selected.features, device.open, "QC")
   print.sigMatrix(cuff, device.open, feature.level)
   print.MDSplot(cuff, selected.features, device.open, show.replicates, log.transform)
   print.PCAplot(selected.features, device.open, "QC", show.replicates)
   print.DistHeat(selected.features, device.open, "QC", show.replicates, log.transform)
   print.dendrogram(selected.features, device.open, "QC", show.replicates, log.transform)

   # Generate plots for all pair-wise sample comparisons
   samples <- samples(cuff@genes)
   count <- NROW(samples)
   if (count == 1) {
      # Skip these if there is only one item.  This should never happen.
      print("Too few samples; skipping pair-wise plots") 
   }
   else {
      for (i in 1:(count-1)) {
         for (j in (i+1):count) {
            currI <- samples[i]
            currJ <- samples[j]
            print.volcanoPlot(selected.features, currI, currJ, device.open, "QC")
            # Plot from next line is identical to the above even with args reversed.  This is an acknowledged bug in cummeRbund. 
            #print.volcanoPlot(selected.features, currJ, currI, device.open, "QC")
            print.scatterPlot(selected.features, currI, currJ, device.open, "QC", log.transform)
            print.scatterPlot(selected.features, currJ, currI, device.open, "QC", log.transform)
            print.MAplot(selected.features, currI, currJ, device.open, "QC", log.transform)
            print.MAplot(selected.features, currJ, currI, device.open, "QC", log.transform)
         }
      }
   }
}

write.significance.data <- function(data, accessor.function, report.name) {
   tryCatch({
      feature_data <- accessor.function(data)
      sig_data <- subset(feature_data, (significant == 'yes'))
      if (nrow(sig_data) > 0) {
         write.table(sig_data, report.name, sep='\t', row.names = F, col.names = T, quote = F)
      }
      else {
         print(paste0("Skipping ", report_name, " - no data found"))
      }
   },
   error = function(err) {
      print(paste0("Error printing ", report.name, " - skipping"))
      print(err)
   })
}

print.dispersionPlot <- build.standardPlotter("Dispersion", 
   function(selected.features, show.replicates, log.transform) {
      return(dispersionPlot(selected.features))
   }
)

print.fpkmSCVPlot <- build.standardPlotter("FPKM.SCV", 
   function(selected.features, show.replicates, log.transform) {
      return(fpkmSCVPlot(selected.features))
   }
)

print.densityPlot <- build.standardPlotter("Density", 
   function(selected.features, show.replicates, log.transform) {
      return(csDensity(selected.features, replicates=show.replicates, logMode=log.transform))
   }
)

print.Boxplot <- build.standardPlotter("Boxplot", 
   function(selected.features, show.replicates, log.transform) {
      return(csBoxplot(selected.features, replicates=show.replicates, logMode=log.transform))
   }
)

print.scatterMatrix <- build.standardPlotter("Boxplot", 
   function(selected.features, show.replicates, log.transform) {
      return(csScatterMatrix(selected.features, replicates=show.replicates, logMode=log.transform))
   }
)

print.volcanoMatrix <- build.standardPlotter("VolcanoMatrix", 
   function(selected.features, show.replicates, log.transform) {
      return(csVolcanoMatrix(selected.features))
   }
)

print.DistHeat <- build.standardPlotter("JSDistanceHeatmap.Samples", 
   function(selected.features, show.replicates, log.transform) {
      return(csDistHeat(selected.features, replicates=show.replicates, samples.not.genes=TRUE, logMode=log.transform))
   }
)

print.PCAplot <- build.standardPlotter("DimensionalityReduction.pca", 
   function(selected.features, show.replicates, log.transform) {
      return(PCAplot(selected.features, replicates=show.replicates))
   }
)

print.sigMatrix <- function(cuff, device.open, feature.level) {
   plotname <- paste0("QC.SignificanceMatrix")
   if (NROW(samples(cuff)) < 3) {
      print(paste0("Skipping the ", plotname, " plot - too few samples"))
   }
   else {
      tryCatch({
         # The significance matrix behaves differently in that it handles the selection within the call.
         sM <- sigMatrix(cuff, level=feature.level)
         print.plotObject(sM, plotname, device.open)
      },
      error = function(err) {
         print(paste0("Error printing the ", plotname, " plot - skipping"))
         print(err)
      })
   }
}

print.MDSplot.unguarded <- build.standardPlotter("DimensionalityReduction.mds", 
   function(selected.features, show.replicates, log.transform) {
      return(MDSplot(selected.features, replicates=show.replicates, logMode=log.transform))
   }
)

print.MDSplot <- function(cuff, selected.features, device.open, show.replicates, log.transform) {
   # The MDSplot throws an error when trying to plot fewer than two samples/replicates.
   if ((!show.replicates && NROW(samples(cuff)) < 3) || (show.replicates && NROW(replicates(cuff)) < 3)) {
      print(paste0("Skipping the QC.DimensionalityReduction.mds plot - too few samples/replicates"))
   }
   else {
      print.MDSplot.unguarded(selected.features, device.open, "QC", show.replicates, log.transform)
   }
}
