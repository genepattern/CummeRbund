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
   print(genesetIds)
   
   feature.selector <- get.feature.selector(feature.level)

   cuff <- readCufflinks.silent(cuffdiff.job, gtf.file, genome.file)
   
   geneset <- getGenes(cuff, genesetIds)
   geneset

   selected.features <- geneset
   if (feature.level != "genes") { selected.features <- feature.selector(geneset) }

   # Write out the various plots
   print.heatmap(selected.features, device.open, "GeneSet", show.replicates, log.transform)
   print.expressonBarplot(selected.features, device.open, "GeneSet", show.replicates, log.transform)
   print.expressonPlot(selected.features, device.open, "GeneSet", show.replicates, log.transform)
   print.dendrogram(selected.features, device.open, "GeneSet", show.replicates, log.transform)
   print.clusterPlot(selected.features, device.open, cluster.count, show.replicates, log.transform)
   
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
            print.volcanoPlot(selected.features, currI, currJ, device.open, "GeneSet")
            # Plot from next line is identical to the above even with args reversed.  This is an acknowledged bug in cummeRbund. 
            #print.volcanoPlot(selected.features, currJ, currI, device.open, "GeneSet")
            print.scatterPlot(selected.features, currI, currJ, device.open, "GeneSet", log.transform)
            print.scatterPlot(selected.features, currJ, currI, device.open, "GeneSet", log.transform)
         }
      }
   }
}

print.heatmap <- build.standardPlotter("Heatmap", 
   function(selected.features, show.replicates, log.transform) {
      return(csHeatmap(selected.features, cluster='both', replicates=show.replicates, logMode=log.transform))
   }
)

print.clusterPlot <- function(selected.features, device.open, cluster.count, show.replicates, log.transform) {
   if (is.null(cluster.count)) {
      print("No cluster.count specified; skipping GeneSet.ClusterPlot") 
   }
   else {
      plotname <- paste0("GeneSet.ClusterPlot.k_",cluster.count)
      tryCatch({
         ic <- csCluster(selected.features, k=cluster.count, logMode=log.transform)
         icp<-csClusterPlot(ic)
         print.plotObject(icp, plotname, device.open)
      },
      error = function(err) {
         print(paste0("Error printing the ", plotname, " plot - skipping"))
         print(err)
      })
   }
}
