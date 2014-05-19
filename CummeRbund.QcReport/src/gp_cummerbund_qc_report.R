## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GP.CummeRbund.QC.Report <- function(cuffdiff.job, gtf.file, genome, output.format,
                                    feature.level, report.as.aggregate, log.transform,
                                    pca.x, pca.y) {
   use.replicates <- !report.as.aggregate
   
   device.open <- get.device.open(output.format)
   
   feature.selector <- get.feature.selector(feature.level)

   cuff <- readCufflinks.silent(cuffdiff.job, gtf.file, genome)
   selected.features <- feature.selector(cuff)
   
   # Write out a table of the differentially expressed features
   write.significance.data(diffData(selected.features), 
                           paste0('QC.sig_diffExp_', feature.level,'.txt'))

   # Write out tables of the significant promoters, splicing, and relCDS data
   write.significance.data(distValues(promoters(cuff)), 'QC.sig_promoter_data.txt')
   write.significance.data(distValues(splicing(cuff)), 'QC.sig_splicing_data.txt')
   write.significance.data(distValues(relCDS(cuff)), 'QC.sig_relCDS_data.txt')

   # Write out the various plots
   print.dispersionPlot(selected.features, device.open, "QC")
   print.fpkmSCVPlot(selected.features, device.open, "QC")
   print.densityPlot(selected.features, device.open, "QC", use.replicates, log.transform)
   print.Boxplot(selected.features, device.open, "QC", use.replicates, log.transform)
   print.MDSplot(cuff, selected.features, device.open, use.replicates, log.transform)
   print.PCAplot(selected.features, pca.x, pca.y, device.open, "QC", use.replicates, log.transform)
   print.DistHeat(selected.features, device.open, "QC", use.replicates, log.transform)
   print.dendrogram(selected.features, device.open, "QC", use.replicates, log.transform)
}

write.significance.data <- function(data, report.name) {
   tryCatch({
      sig_data <- subset(data, (significant == 'yes'))
      if (nrow(sig_data) > 0) {
         write.table(sig_data, report.name, sep='\t', row.names = F, col.names = T, quote = F)
      }
      else {
         print(paste0("Skipping ", report.name, " - no data found meeting significance threshold"))
      }
   },
   error = function(err) {
      print(paste0("Error printing ", report.name, " - skipping"))
      print(conditionMessage(err))
   })
}

print.dispersionPlot <- build.standardPlotter("Dispersion", 
   function(selected.features, use.replicates, log.transform) {
      return(dispersionPlot(selected.features))
   }
)

print.fpkmSCVPlot <- build.standardPlotter("FPKM.SCV", 
   function(selected.features, use.replicates, log.transform) {
      return(fpkmSCVPlot(selected.features))
   }
)

print.densityPlot <- build.standardPlotter("Density", 
   function(selected.features, use.replicates, log.transform) {
      return(csDensity(selected.features, replicates=use.replicates, logMode=log.transform))
   }
)

print.Boxplot <- build.standardPlotter("Boxplot", 
   function(selected.features, use.replicates, log.transform) {
      return(csBoxplot(selected.features, replicates=use.replicates, logMode=log.transform))
   }
)

print.DistHeat <- build.standardPlotter("JSDistanceHeatmap.Samples", 
   function(selected.features, use.replicates, log.transform) {
      return(csDistHeat(selected.features, replicates=use.replicates, samples.not.genes=TRUE, logMode=log.transform))
   }
)

print.PCAplot <- build.XYAxisPlotter("DimensionalityReduction.pca", 
   function(selected.features, x, y, use.replicates, log.transform) {
      return(PCAplot(selected.features, x, y, replicates=use.replicates))
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
         print(conditionMessage(err))
      })
   }
}

print.MDSplot.unguarded <- build.standardPlotter("DimensionalityReduction.mds", 
   function(selected.features, use.replicates, log.transform) {
      return(MDSplot(selected.features, replicates=use.replicates, logMode=log.transform))
   }
)

print.MDSplot <- function(cuff, selected.features, device.open, use.replicates, log.transform) {
   # The MDSplot throws an error when trying to plot fewer than two samples/replicates.
   if ((!use.replicates && NROW(samples(cuff)) < 3) || (use.replicates && NROW(replicates(cuff)) < 3)) {
      print(paste0("Skipping the QC.DimensionalityReduction.mds plot - too few samples/replicates"))
   }
   else {
      print.MDSplot.unguarded(selected.features, device.open, "QC", use.replicates, log.transform)
   }
}

check.pca.axis.selection <- function(axis_selection) {
   if (!grepl('^PC\\d*$', axis_selection)) {
      stop(paste0("Unrecognized PCA component axis selector '", axis_selection, "'.  Must be of the form PC1, PC2, PC3, etc."))
   }
}
