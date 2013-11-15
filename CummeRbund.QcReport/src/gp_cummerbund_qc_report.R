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
   #ggsave(build.filename("DispersionPlot", output.format), plot=disp)
   
   genes.scv<-fpkmSCVPlot(genes(cuff))
   print.plotObject(genes.scv, "Genes.fpkmSCVPlot", device.open)
   #ggsave(build.filename("Genes.fpkmSCVPlot", output.format),plot=genes.scv)
   
   isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
   print.plotObject(isoforms.scv, "Isoforms.fpkmSCVPlot", device.open) 
   #ggsave(build.filename("Isoforms.fpkmSCVPlot", output.format), plot=isoforms.scv)
   
   dens<-csDensity(genes(cuff))
   print.plotObject(dens, "Density", device.open) 
   #ggsave(build.filename("Density", output.format), plot=dens)

   densRep<-csDensity(genes(cuff),replicates=T)
   print.plotObject(densRep, "Density.rep", device.open) 
   #ggsave(build.filename("Density.rep", output.format), plot=densRep)

   b<-csBoxplot(genes(cuff))
   print.plotObject(b, "Boxplot", device.open) 
   #ggsave(build.filename("Boxplot", output.format), plot=b)

   brep<-csBoxplot(genes(cuff),replicates=T)
   print.plotObject(brep, "Boxplot.rep", device.open) 
   #ggsave(build.filename("Boxplot.rep", output.format), plot=brep)

   s<-csScatterMatrix(genes(cuff))
   print.plotObject(s, "ScatterMatrix", device.open) 
   #ggsave(build.filename("ScatterMatrix", output.format), plot=s)

   v<-csVolcanoMatrix(genes(cuff))
   print.plotObject(v, "VolcanoMatrix", device.open) 
   #ggsave(build.filename("VolcanoMatrix", output.format), plot=v)

   # Not yet: need to figure out how to pass labels
   # s<-csScatter(genes(cuff),"hESC","Fibroblasts",smooth=T)
   # m<-MAplot(genes(cuff),"hESC","Fibroblasts")
   # mCount<-MAplot(genes(cuff),"hESC","Fibroblasts",useCount=T)

   # The dendrogram plots behave differently - apparently the print/plot call is embedded within.
   device.open("Dendrogram")
   #pdf(build.filename("Dendrogram", "pdf"))
   #dend<-csDendro(genes(cuff))
   dev.off()
   #ggsave(build.filename("Dendrogram", output.format))

   device.open("Dendrogram.rep")
   #pdf(build.filename("Dendrogram.rep", "pdf"))
   #dendRep<-csDendro(genes(cuff),replicates=T)
   dev.off()
   #ggsave(build.filename("Dendrogram.rep", output.format))
}

get.device.open <- function(extension) {
   if (extension == "pdf") { 
      return(function(filename_base) {
         pdf(paste0(filename_base, ".pdf"))
      })
   }
   if (extension == "svg") { 
      return(function(filename_base) {
         svg(paste0(filename_base, ".svg"))
      })
   }
   if (extension == "png") { 
      return(function(filename_base) {
         png(paste0(filename_base, ".png"))
      })
   }
   stop("Unhandled extension")
}

print.plotObject <- function(plotObj, filename_base, device.open) {
   device.open(filename_base)
   print(plotObj)
   dev.off()
}