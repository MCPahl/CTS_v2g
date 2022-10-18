library(tidyverse)

files = list.files()
files = files[grepl("\\.results", files)]
files = files[grepl("CTS", files)]

ldsr = lapply(files, function(f){
	dat <- read.delim(f) %>%
	filter(row_number()==1) %>% 
	mutate(cell = gsub("_promoterContactingOCRs.*", "", f))
	dat
})

ldsr = do.call("rbind", ldsr)
pldsr = ldsr %>% 
	mutate(p.adjust = p.adjust(ldsr$Enrichment_p, method = "bonferroni")) %>%
	mutate(sig = ifelse(p.adjust < 0.05, "sig", ifelse(Enrichment_p < 0.05, "nomimal", "ns"))) %>% 
    mutate(Category = "Carpal Tunnel Syndrome") %>% 
    mutate(z = -qnorm(Enrichment_p))


pdf("pldsr_output.pdf", useDingbats=FALSE, height=5, width=5)
 ggplot( pldsr, aes(y= cell, x = Category))+
    geom_point(aes(size= Enrichment, col= z))+
    scale_color_viridis(discrete=FALSE)+
    geom_text(aes(label = sig), size = 2, color = "white") +
    theme_bw()
dev.off()




pdf("pldsr_enrichment_results.pdf", useDingbats=FALSE)
ggplot(pldsr, aes(x=cell, y=Enrichment, col=sig, fill=sig))+
 geom_point(size=5)+
 geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error*1.96, ymax=Enrichment+Enrichment_std_error*1.96), width=0)+
 scale_colour_manual(values = c("orange","black", "red"))+
 geom_hline(yintercept = 1, linetype="dotted", color="blue", size=1)+
 coord_flip()+
 theme_bw()+
  scale_x_discrete(limits = rev(levels(pldsr$Category)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
 dev.off()

write.table(pldsr, file = "pldsr_results.txt", quote=F, sep="\t",row.names=F)