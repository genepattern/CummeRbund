## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GP.CummeRbund.QC.Report <- function(cuffdiff.job, gtf.file, genome.file, output.format) {
   library(cummeRbund)

   device.open <- get.device.open(output.format)

   ## Calls based on examples from http://compbio.mit.edu/cummeRbund/manual_2_0.html
   cuff<-readCufflinks(cuffdiff.job, gtfFile = gtf.file, genome = genome.file)
   cuff
   
   disp<-dispersionPlot(genes(cuff))
   print.plotObject(disp, "DispersionPlot", device.open)
   
   genes.scv<-fpkmSCVPlot(genes(cuff))
   print.plotObject(genes.scv, "Genes.fpkmSCVPlot", device.open)
   
   isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
   print.plotObject(isoforms.scv, "Isoforms.fpkmSCVPlot", device.open) 
   
   dens<-csDensity(genes(cuff))
   print.plotObject(dens, "Density", device.open) 

   densRep<-csDensity(genes(cuff),replicates=T)
   print.plotObject(densRep, "Density.rep", device.open) 

   b<-csBoxplot(genes(cuff))
   print.plotObject(b, "Boxplot", device.open) 

   brep<-csBoxplot(genes(cuff),replicates=T)
   print.plotObject(brep, "Boxplot.rep", device.open) 

   s<-csScatterMatrix(genes(cuff))
   print.plotObject(s, "ScatterMatrix", device.open) 

   v<-csVolcanoMatrix(genes(cuff))
   print.plotObject(v, "VolcanoMatrix", device.open) 

   # Not yet: need to figure out how to pass labels
   # s<-csScatter(genes(cuff),"hESC","Fibroblasts",smooth=T)
   # m<-MAplot(genes(cuff),"hESC","Fibroblasts")
   # mCount<-MAplot(genes(cuff),"hESC","Fibroblasts",useCount=T)

   # The dendrogram plots behave differently - apparently the print/plot call is embedded within.
   device.open("Dendrogram")
   dend<-csDendro(genes(cuff))
   dev.off()

   device.open("Dendrogram.rep")
   dendRep<-csDendro(genes(cuff),replicates=T)
   dev.off()
}
