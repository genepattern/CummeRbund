## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GP.CummeRbund.SelectedGene.Report <- function(cuffdiff.job, feature.id, gtf.file, genome.file,
                                              output.format, feature.level, show.replicates, log.transform) {
   device.open <- get.device.open(output.format)
   
   feature.selector <- get.feature.selector(feature.level)

   cuff <- readCufflinks.silent(cuffdiff.job, gtf.file, genome.file)
   checkCuffVersionAbove2(cuff)
   print(cuff)

   print(paste0("Looking up the selected gene using feature ID '", feature.id, "'"))
   selected.gene <- getGene(cuff, feature.id)
   print(selected.gene)

   selected.features <- selected.gene
   if (feature.level != "genes") { selected.features <- feature.selector(selected.gene) }

   print.expressonBarplot(selected.features, show.replicates, log.transform, device.open, "SelectedGene")
   print.expressonPlot(selected.features, show.replicates, log.transform, device.open, "SelectedGene")

   tryCatch({
      mySimilar <- findSimilar(cuff, selected.gene, n=5)
      mySimilar.expression<-expressionPlot(mySimilar, logMode=log.transformm, showErrorbars=F)
      print.plotObject(mySimilar.expression, "SelectedGene.SimilarityExpressionPlot", device.open)
   },
   error = function(err) {
      print("Error printing the SelectedGene.SimilarityExpressionPlot plot - skipping")
      print(err)
   })
}
