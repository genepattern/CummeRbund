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
   # Verify we have enough conditions for plotting
   conditions.count <- NROW(selected.conditions)
   if (conditions.count < 2) {
      stop("Too few conditions for plotting.  Must have two or more.") 
   }

   use.replicates <- !report.as.aggregate
   device.open <- get.device.open(output.format)

   feature.selector <- get.feature.selector(feature.level)

   cuff <- readCufflinks.silent(cuffdiff.job, gtf.file, genome)
   
   check.selected.conditions(selected.conditions, cuff)
   selected.data <- feature.selector(cuff)
      
   # Generate plots for pairwise sample comparisons
   for (i in 1:(conditions.count-1)) {
      for (j in (i+1):conditions.count) {
         currI <- selected.conditions[i]
         currJ <- selected.conditions[j]
         print.volcanoPlot(selected.data, currI, currJ, device.open, "SelectedConditions")
         # Plot from next line is identical to the above even with args reversed.  This is an acknowledged bug in cummeRbund. 
         #print.volcanoPlot(selected.data, currJ, currI, device.open, "SelectedConditions")
         print.scatterPlot(selected.data, currI, currJ, device.open, "SelectedConditions", log.transform)
         print.scatterPlot(selected.data, currJ, currI, device.open, "SelectedConditions", log.transform)
         print.MAplot(selected.data, currI, currJ, device.open, "SelectedConditions", log.transform)
         print.MAplot(selected.data, currJ, currI, device.open, "SelectedConditions", log.transform)
      }
   }
}
