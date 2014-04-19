## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GP.CummeRbund.SelectedConditions.Report <- function(cuffdiff.job, selected.conditions, gtf.file, genome, output.format, 
                                                    feature.level, report.as.aggregate, log.transform) {
   use.replicates <- !report.as.aggregate
   device.open <- get.device.open(output.format)

   feature.selector <- get.feature.selector(feature.level)

   cuff <- readCufflinks.silent(cuffdiff.job, gtf.file, genome)
   
   check.selected.conditions(selected.conditions, cuff)

   # Get a CuffGeneSet based on the selected conditions and *all* genes (features).
   #
   # Possible future improvement: modify the underlying queries to work without the geneIdList.
   #   That is quite complex but it would be much more efficient.  If the following works for
   #   fairly large datasets, however, it's probably not worth the effort to make this change
   #   as it will likely become a moot point when we re-implement the workflow with the new
   #   Cufflinks 2.2.x improvements for working with subsets of conditions upstream.    
   genesetForConditions <- getGenes(cuff, geneIdList=featureNames(cuff@genes), sampleIdList=selected.conditions)

   selected.data <- genesetForConditions
   if (feature.level != "genes") { selected.data <- feature.selector(genesetForConditions) }

   # Write out the various plots
   print.heatmap(selected.data, device.open, "SelectedConditions", use.replicates, log.transform)
   print.expressonBarplot(selected.data, device.open, "SelectedConditions", use.replicates, log.transform)
   print.expressonPlot(selected.data, device.open, "SelectedConditions", use.replicates, log.transform)
   print.dendrogram(selected.data, device.open, "SelectedConditions", use.replicates, log.transform)
   
   # Generate plots for pairwise sample comparisons
   # This is skipped if no conditions were specified in selected.conditions due to the overhead and
   # potential large number of files generated.  
   if (conditions.count < 1) {
      # Skip these if there is not at least one pair
      print("Too few conditions; skipping pairwise plots") 
   }
   else {
      for (i in 1:(conditions.count-1)) {
         for (j in (i+1):conditions.count) {
            currI <- selected.conditions[i]
            currJ <- selected.conditions[j]
            print.volcanoPlot(selected.data, currI, currJ, device.open, "SelectedConditions")
            # Plot from next line is identical to the above even with args reversed.  This is an acknowledged bug in cummeRbund. 
            #print.volcanoPlot(selected.data, currJ, currI, device.open, "SelectedConditions")
            print.scatterPlot(selected.data, currI, currJ, device.open, "SelectedConditions", log.transform)
            print.scatterPlot(selected.data, currJ, currI, device.open, "SelectedConditions", log.transform)
         }
      }
   }
}

# Move to Util
print.heatmap <- build.standardPlotter("Heatmap", 
   function(selected.data, use.replicates, log.transform) {
      return(csHeatmap(selected.data, clustering='both', replicates=use.replicates, logMode=log.transform))
   }
)