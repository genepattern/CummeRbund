## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GP.CummeRbund.Geneset.Report <- function(cuffdiff.job, geneset.file, selected.conditions, gtf.file, genome,
                                         output.format, feature.level, report.as.aggregate, log.transform) {
   use.replicates <- !report.as.aggregate
   device.open <- get.device.open(output.format)

   genesetIds <- scan(geneset.file, what="character", sep="\n", quiet=TRUE)
   
   print("Generating plots based on the following genes:")
   print(genesetIds)

   feature.selector <- get.feature.selector(feature.level)

   cuff <- readCufflinks.silent(cuffdiff.job, gtf.file, genome)
   
   conditions.count <- NROW(selected.conditions)
   if (conditions.count == 0) {
      selected.conditions <- NULL
   }
   else {
      check.selected.conditions(selected.conditions, cuff)
   }
   
   geneset <- getGenes(cuff, genesetIds, sampleIdList=selected.conditions)
   geneset

   selected.features <- geneset
   if (feature.level != "genes") { selected.features <- feature.selector(geneset) }

   # Write out the various plots
   print.heatmap(selected.features, device.open, "GeneSet", use.replicates, log.transform)
   print.expressonBarplot(selected.features, device.open, "GeneSet", use.replicates, log.transform)
   print.expressonPlot(selected.features, device.open, "GeneSet", use.replicates, log.transform)
   print.dendrogram(selected.features, device.open, "GeneSet", use.replicates, log.transform)
   
   # Generate plots for pairwise sample comparisons
   # This is skipped if no conditions were specified in selected.conditions due to the overhead and
   # potential large number of files generated.  
   if (conditions.count < 2) {
      # Skip these if there is not at least one pair
      print("Too few conditions; skipping pairwise plots") 
   }
   else {
      for (i in 1:(conditions.count-1)) {
         for (j in (i+1):conditions.count) {
            currI <- selected.conditions[i]
            currJ <- selected.conditions[j]
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
   function(selected.features, use.replicates, log.transform) {
      return(csHeatmap(selected.features, clustering='both', replicates=use.replicates, logMode=log.transform))
   }
)