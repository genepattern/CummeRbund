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
   print(cuff)
   features <- feature.selector(cuff)
   
   # Write out a table of the differentially expressed features
   feature_diff_data <- diffData(features)
   sig_feature_data <- subset(feature_diff_data, (significant == 'yes'))
   write.table(sig_feature_data, paste0('sig_diffExp_', feature.level,'.txt'), sep='\t',
               row.names = F, col.names = T, quote = F)

   ## Calls based on examples from http://compbio.mit.edu/cummeRbund/manual_2_0.html
   disp <- dispersionPlot(features)
   print.plotObject(disp, "Dispersion", device.open)
   
   fpkm.scv <- fpkmSCVPlot(features)
   print.plotObject(fpkm.scv, "FPKM.SCV", device.open)

   dens <- csDensity(features, replicates=show.replicates, logMode=log.transform)
   print.plotObject(dens, "Density", device.open) 

   box <- csBoxplot(features, replicates=show.replicates, logMode=log.transform)
   print.plotObject(box, "Boxplot", device.open) 

   s <- csScatterMatrix(features, replicates=show.replicates, logMode=log.transform)
   print.plotObject(s, "Scatter", device.open) 

   v <- csVolcanoMatrix(features)
   print.plotObject(v, "Volcano", device.open) 

   # The significance matrix behaves differently in that it handles the selection within the call.
   sM <- sigMatrix(cuff, level=feature.level)
   print.plotObject(sM, "SignificanceMatrix", device.open) 

   # The MDSplot throws an error when trying to plot fewer than two samples/replicates.
   if ((!show.replicates && NROW(samples(cuff)) < 3) || (show.replicates && NROW(replicates(cuff)) < 3)) {
      print("Skipping the MDS plot - too few samples/replicates")
   }
   else {
      mds <- MDSplot(features, replicates=show.replicates, logMode=log.transform)
      print.plotObject(mds, "DimensionalityReduction.mds", device.open)
   }

   pca <- PCAplot(features, replicates=show.replicates)
   print.plotObject(pca, "DimensionalityReduction.pca", device.open) 

   distHeatSamples <- csDistHeat(features, replicates=show.replicates, samples.not.genes=TRUE, logMode=log.transform)
   print.plotObject(distHeatSamples, "JSDistanceHeatmap.Samples", device.open) 

   # Not yet: need to figure out how to pass labels
   # s<-csScatter(features,"hESC","Fibroblasts",smooth=T)
   # m<-MAplot(features,"hESC","Fibroblasts")
   # mCount<-MAplot(features,"hESC","Fibroblasts",useCount=T)

   # The dendrogram plots behave differently - apparently the print/plot call is embedded within.
   device.open("Dendrogram")
   dend <- csDendro(features,replicates=show.replicates,logMode=log.transform)
   dev.off()
}
