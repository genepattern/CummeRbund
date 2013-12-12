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
   checkCuffVersionAbove2(cuff)
   print(cuff)
   selected.features <- feature.selector(cuff)
   
   # Write out a table of the differentially expressed features
   tryCatch({
      feature_diff_data <- diffData(selected.features)
      sig_feature_data <- subset(feature_diff_data, (significant == 'yes'))
      if (nrow(sig_feature_data) > 0) {
         write.table(sig_feature_data, paste0('QC.sig_diffExp_', feature.level,'.txt'), sep='\t',
                     row.names = F, col.names = T, quote = F)
      }
      else {
         print("Skipping differentially expressed feature table - no data found")
      }
   },
   error = function(err) {
      print("Error printing table of the differentially expressed features - skipping")
      print(err)
   })

   # Write out tables of the significant promoters, splicing, and relCDS data
   tryCatch({
      promoter_diff_data <- distValues(promoters(cuff))
      sig_promoter_data <- subset(promoter_diff_data, (significant == 'yes')) 
      if (nrow(sig_promoter_data) > 0) {
         write.table(sig_promoter_data, 'QC.sig_promoter_data.txt', sep='\t',
                     row.names = F, col.names = T, quote = F)
      }
      else {
         print("Skipping significant promoter table - no data found")
      }
   },
   error = function(err) {
      print("Error printing table of the significant promoter data - skipping")
      print(err)
   })
   tryCatch({
      splicing_diff_data <- distValues(splicing(cuff))
      sig_splicing_data <- subset(splicing_diff_data, (significant == 'yes'))
      if (nrow(sig_promoter_data) > 0) {
         write.table(sig_splicing_data, 'QC.sig_splicing_data.txt', sep='\t',
                     row.names = F, col.names = T, quote = F)
      }
      else {
         print("Skipping significant splicing table - no data found")
      }
   },
   error = function(err) {
      print("Error printing table of the significant splicing data - skipping")
      print(err)
   })
   tryCatch({
      relCDS_diff_data <- distValues(relCDS(cuff))
      sig_relCDS_data <- subset(relCDS_diff_data, (significant == 'yes'))
      if (nrow(sig_relCDS_data) > 0) {
         write.table(sig_relCDS_data, 'QC.sig_relCDS_data.txt', sep='\t',
                     row.names = F, col.names = T, quote = F)
      }
      else {
         print("Skipping significant relCDS table - no data found")
      }
   },
   error = function(err) {
      print("Error printing table of the significant relCDS data - skipping")
      print(err)
   })
   
   tryCatch({
      disp <- dispersionPlot(selected.features)
      print.plotObject(disp, "QC.Dispersion", device.open)
   },
   error = function(err) {
      print("Error printing the QC.Dispersion plot - skipping")
      print(err)
   })
   
   tryCatch({
      fpkm.scv <- fpkmSCVPlot(selected.features)
      print.plotObject(fpkm.scv, "QC.FPKM.SCV", device.open)
   },
   error = function(err) {
      print("Error printing the QC.FPKM.SCV plot - skipping")
      print(err)
   })

   tryCatch({
      dens <- csDensity(selected.features, replicates=show.replicates, logMode=log.transform)
      print.plotObject(dens, "QC.Density", device.open) 
   },
   error = function(err) {
      print("Error printing the QC.Density plot - skipping")
      print(err)
   })

   tryCatch({
      box <- csBoxplot(selected.features, replicates=show.replicates, logMode=log.transform)
      print.plotObject(box, "QC.Boxplot", device.open) 
   },
   error = function(err) {
      print("Error printing the QC.Boxplot plot - skipping")
      print(err)
   })

   tryCatch({
      s <- csScatterMatrix(selected.features, replicates=show.replicates, logMode=log.transform)
      print.plotObject(s, "QC.ScatterMatrix", device.open) 
   },
   error = function(err) {
      print("Error printing the QC.ScatterMatrix plot - skipping")
      print(err)
   })

   tryCatch({
      v <- csVolcanoMatrix(selected.features)
      print.plotObject(v, "QC.VolcanoMatrix", device.open) 
   },
   error = function(err) {
      print("Error printing the QC.VolcanoMatrix plot - skipping")
      print(err)
   })

   # The sigMatrix is not useful with fewer than three samples
   if (NROW(samples(cuff)) < 3) {
      print("Skipping the QC.SignificanceMatrix plot - too few samples")
   }
   else {
      tryCatch({
         # The significance matrix behaves differently in that it handles the selection within the call.
         sM <- sigMatrix(cuff, level=feature.level)
         print.plotObject(sM, "QC.SignificanceMatrix", device.open)
      },
      error = function(err) {
         print("Error printing the QC.SignificanceMatrix plot - skipping")
         print(err)
      })
   }

   # The MDSplot throws an error when trying to plot fewer than two samples/replicates.
   if ((!show.replicates && NROW(samples(cuff)) < 3) || (show.replicates && NROW(replicates(cuff)) < 3)) {
      print("Skipping the MDS plot - too few samples/replicates")
   }
   else {
      tryCatch({
         mds <- MDSplot(selected.features, replicates=show.replicates, logMode=log.transform)
         print.plotObject(mds, "QC.DimensionalityReduction.mds", device.open)
      },
      error = function(err) {
         print("Error printing the QC.DimensionalityReduction.mds plot - skipping")
         print(err)
      })
   }

   tryCatch({
      pca <- PCAplot(selected.features, replicates=show.replicates)
      print.plotObject(pca, "QC.DimensionalityReduction.pca", device.open) 
   },
   error = function(err) {
      print("Error printing the QC.DimensionalityReduction.pca plot - skipping")
      print(err)
   })

   tryCatch({
      distHeatSamples <- csDistHeat(selected.features, replicates=show.replicates, samples.not.genes=TRUE, logMode=log.transform)
      print.plotObject(distHeatSamples, "QC.JSDistanceHeatmap.Samples", device.open) 
   },
   error = function(err) {
      print("Error printing the QC.JSDistanceHeatmap.Samples plot - skipping")
      print(err)
   })

   print.dendrogram(selected.features, show.replicates, log.transform, device.open, "QC")

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
            # Seems to be identical to the above even with args reversed.
            #print.volcanoPlot(selected.features, currJ, currI, device.open, "QC")
            print.scatterPlot(selected.features, currI, currJ, log.transform, device.open, "QC")
            print.scatterPlot(selected.features, currJ, currI, log.transform, device.open, "QC")
            print.MAplot(selected.features, currI, currJ, log.transform, device.open, "QC")
            print.MAplot(selected.features, currJ, currI, log.transform, device.open, "QC")
         }
      }
   }
}
