## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GP.CummeRbund.SelectedGene.Report <- function(cuffdiff.job, feature.id, selected.conditions, find.similar, gtf.file,
                                              genome, output.format, feature.level, report.as.aggregate, log.transform) {
   use.replicates <- !report.as.aggregate
   device.open <- get.device.open(output.format)
   
   feature.selector <- get.feature.selector(feature.level)

   cuff <- readCufflinks.silent(cuffdiff.job, gtf.file, genome)

   conditions.count <- NROW(selected.conditions)
   if (conditions.count == 0) {
      selected.conditions <- NULL
   }
   else {
      check.selected.conditions(selected.conditions, cuff)
   }

   print(paste0("Looking up the selected gene using feature ID '", feature.id, "'"))
   selected.gene <- getGene(cuff, feature.id, sampleIdList=selected.conditions)
   print(selected.gene)

   selected.features <- selected.gene
   if (feature.level != "genes") { selected.features <- feature.selector(selected.gene) }

   print.expressonBarplot(selected.features, device.open, "SelectedGene", use.replicates, log.transform)
   print.expressonPlot(selected.features, device.open, "SelectedGene", use.replicates, log.transform)

   if (is.null(find.similar)) {
      print("find.similar not specified; skipping SelectedGene.SimilarityExpressionPlot")
   }
   else {
      similar.features <- findSimilar(cuff, feature.id, n=find.similar)
      print.similarityPlot(similar.features, device.open, "SelectedGene", log.transform=log.transform)
      print.similarityBarplot(similar.features, device.open, "SelectedGene", log.transform=log.transform)
   }
}

print.similarityPlot <- build.standardPlotter("SimilarityExpressionPlot", 
   function(selected.features, use.replicates, log.transform) {
      return(expressionPlot(selected.features, logMode=log.transform, showErrorbars=FALSE))
   }
)

print.similarityBarplot <- build.standardPlotter("SimilarityExpressionBarplot", 
   function(selected.features, use.replicates, log.transform) {
      return(expressionBarplot(selected.features, logMode=log.transform, showErrorbars=FALSE))
   }
)
