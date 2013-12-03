## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GP.CummeRbund.Geneset.Report <- function(cuffdiff.job, geneset.file, gtf.file, genome.file, output.format,
                                         feature.level, show.replicates, log.transform, cluster.count) {
   device.open <- get.device.open(output.format)

   genesetIds <- unlist(read.table(geneset.file))
   
   print("Generating plots based on the following genes:")
   genesetIds
   
   feature.selector <- get.feature.selector(feature.level)

   cuff <- readCufflinks.silent(cuffdiff.job, gtf.file, genome.file)
   cuff
   
   geneset <- getGenes(cuff, genesetIds)
   geneset

   features <- geneset
   if (feature.level != "genes") { features <- feature.selector(geneset) }

   ## Calls based on examples from http://compbio.mit.edu/cummeRbund/manual_2_0.html
   heatmap <- csHeatmap(features, cluster='both', replicates=show.replicates, logMode=log.transform)
   print.plotObject(heatmap, "Heatmap", device.open)

   expBarplot<-expressionBarplot(features, replicates=show.replicates, logMode=log.transform)
   print.plotObject(expBarplot, "ExpressionBarplot", device.open)
   
   expPlot<-expressionPlot(features, replicates=show.replicates, logMode=log.transform)
   print.plotObject(expPlot, "ExpressionPlot", device.open)

   # The dendrogram plots behave differently - apparently the print/plot call is embedded within.
   device.open("Dendrogram")
   dend <- csDendro(features, replicates=show.replicates, logMode=log.transform)
   dev.off()
   
   if (!is.null(cluster.count)) {
      ic <- csCluster(features, k=cluster.count, logMode=log.transform)
      icp<-csClusterPlot(ic)
      print.plotObject(icp, paste0("ClusterPlot.k_",cluster.count), device.open)
   }
   
   # Generate plots for all pair-wise sample comparisons
   samples <- samples(cuff@genes)
   count <- NROW(samples)
   if (count == 1) {
      # Bail out here if there is only one item.  This should never happen
      print("Too few samples; skipping pair-wise plots") 
      return() 
   }
   
   for (i in 1:(count-1)) {
      for (j in (i+1):count) {
         currI <- samples[i]
         currJ <- samples[j]
         v1<-csVolcano(features, x=currI, y=currJ, showSignificant=TRUE)
         print.plotObject(v1, paste0("Volcano.",currI,"_",currJ), device.open)
         v2<-csVolcano(features, x=currJ, y=currI, showSignificant=TRUE)
         print.plotObject(v2, paste0("Volcano.",currJ,"_",currI), device.open)
         s1 <- csScatter(features, x=currI, y=currJ, colorByStatus=TRUE, logMode=log.transform)
         print.plotObject(s1, paste0("Scatter.",currI,"_",currJ), device.open)
         s2 <- csScatter(features, x=currJ, y=currI, colorByStatus=TRUE, logMode=log.transform)
         print.plotObject(s2, paste0("Scatter.",currJ,"_",currI), device.open)
      }
   }
}
