
library(qqman)
library(data.table)
library(tidyverse)

gwas = fread("CTS2022_summarystats/summarystats_CTS_combined_EAF_estimated.txt") %>%
	rename(CHR = Chr, BP = Position_hg38, SNP = rsID) %>%
	mutate(CHR = as.integer(gsub("chr", "", CHR)), SNP = gsub(",.*", "", SNP))

leadSNPs = read.delim("Skuladottir2022_carpaltunnel_leadSNPs.txt")
sentinels = leadSNPs$rs_name


v2g = read.delim("v2g_carpal_tunnel_sfg_database.txt")
muscle.sentinel = unique(v2g[v2g$cell == 'PHMB_Myotube',]$sentinel)


pdf("plots/manhattan_annotated_leadSNPs.pdf")
manhattan(gwas, highlight = sentinels)
dev.off()

png("plots/manhattan_annotated_leadSNPs.png")
manhattan(gwas, highlight = muscle.sentinel)
dev.off()

