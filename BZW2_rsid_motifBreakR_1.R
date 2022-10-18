R/3.6.2

library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MotifDb)
library(TFBSTools)
library(tidyverse)
library(Gviz)
library(motifStack)
proxylist = 'rs2240855'

#### extract variants from dbSNP144
variants <- snps.from.rsid(rsid = proxylist,
                           dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37,
                           search.genome = BSgenome.Hsapiens.UCSC.hg19)

#### subset total MotifDb
hoco=subset(
        MotifDb, dataSource %in% c("HOCOMOCOv10")
        )
hoco=subset(
        hoco, organism %in% c("Hsapiens")
)

#### run motifbreakR
motifbreakr.results <- motifbreakR(snpList = variants, pwmList = hoco, 
                                   filterp = TRUE,
                                   threshold = 1e-3,
                                   method = "ic",
                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25))
 motifbreakr.results = calculatePvalue(motifbreakr.results)

save(motifbreakr.results, 
     file="motifbreakR.RData")

write.table(motifbreakr.results, 
            file="/mnt/isilon/sfgi/pahlm/projects/small_jobs/carpal_tunnel_v2g/motifbreakR.txt",
            sep="\t", row.names=F, quote=F)


motifbreakr.results
pdf("rs197393_BZW2_motif.pdf")
plotMB(results = motifbreakr.results, rsid = "rs197393")
dev.off()

seqnames(motifbreakr.results)[1] = 7


simplified= unique(data.frame(proxy=motifbreakr_result_df$SNP_id, ref=motifbreakr_result_df$REF, alt=motifbreakr_result_df$ALT, gene = motifbreakr_result_df$geneSymbol, diff=motifbreakr_result_df$alleleDiff,  effect =motifbreakr_result_df$alleleDiff))

left_join(dat , simplified)

write.table(snp_in_tfs, 
            file="Desktop/COVID19_SNPs_inTFBS_anno.txt",
            sep="\t", row.names=F, quote=F)


proxyList = c("rs3752495", "rs9932282", "rs8062685", "rs3184470")



"chr16:705360-705360"
names(frag_int)[5:7]=c("chr","start","end")
 x=GRanges(as.data.frame(frag_int))



 pdf("rs197393_BZW2_motif.pdf", height = 20, width=20)
plotMB(results = motifbreakr.results, rsid = "rs2240855")
dev.off()



    motif.starts <- sapply(results$motifPos, `[`, 1)
    motif.starts <- start(results) + motif.starts
    motif.starts <- order(motif.starts)
    results <- results[motif.starts]
    g <- genome(results)[[1]]
    result <- results[names(results) %in% rsid]
    result <- result[order(sapply(result$motifPos, min), sapply(result$motifPos, 
        max)), ]
    result <- result[result$effect %in% effect]
    chromosome <- as.character(seqnames(result))[[1]]
