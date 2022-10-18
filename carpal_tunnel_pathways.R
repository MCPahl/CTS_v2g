library(tidyverse)
library(ggrepel)

v2g =  read.delim("v2g_carpal_tunnel_sfg_database.txt")
muscle_v2g = v2g %>% filter(cell == "PHMB_Myotube") 
osteoblast_v2g = v2g %>% filter(cell == "hMSC_BMP2") 


HG.test = function(grp,anno){
  universe = unique(unlist(anno)) #Set of genes that have been annotated, number of balls 
  gset = grp
  gset = gset[gset %in% universe]
  tmp = lapply(anno,function(j){
      p = length(j[j %in% gset]) #white balls drawn from an urn
      m = length(j) #white balls in the urn
      n = length(universe) - m #Number of black balls in the urn
      k = length(gset) #Number of balls drawn from the urn
      expected = ceiling(length(j)/length(universe)*length(gset)) #Number of white balls expected to be drawn, rounding up (replace ceiling with floor to round down, or remove to leave the decimal)
      hg.pval = phyper(p,m,n,k,lower.tail=FALSE) #Calculate p value
      genes = paste(j[j %in% gset], collapse=":") #Combine the list of genes
      data.frame(p.val = hg.pval,genes= genes, observed=p, expected=expected) #Put the data together in a data.frame
    })
   out = do.call("rbind",tmp) #Merge the list of dataframes into a single dataframe
   out$Pathway = names(tmp) #Convert row.names to a column
   row.names(out) = NULL #Remove the row.names
   out$p.adjust = p.adjust(out$p.val,method="BH") #Calculate p.value after FDR correction
   #out = out[out$p.adjust<0.05 & out$observed>0,]
   out = out[order(out$p.val,decreasing=FALSE), c(5,3,4,1,6,2)]
   out$enrichment = out$observed/out$expected
   out
}

#Reactome
load("/mnt/isilon/sfgi/pahlm/annotationFiles/msigdb_v7.0_GMTs/c5.bp.v7.0.symbols.Rdata")

skeletal_muscle = HG.test(unique(muscle_v2g$gene_name), c5.bp)
osteoblast = HG.test(unique(osteoblast_v2g$gene_name), c5.bp)

reference = c(skeletal_muscle[skeletal_muscle$p.adjust < 0.05,]$Pathway)


pdf("plots/skeletal_muscle_pathways.pdf", useDingbats=FALSE, height = 12, width = 12)
ggplot(skeletal_muscle, aes(x = enrichment, y = -log10(p.adjust)))+
geom_point(aes(size=observed), color = "red")+
theme_bw()+
ylim(0,4)+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
	text = element_text(size = 20))+
	geom_label_repel(data = subset(skeletal_muscle, rank(skeletal_muscle$p.adjust) <= 15), aes(label = Pathway), size =5)
dev.off()


pdf("plots/osteoblast_pathways.pdf", useDingbats=FALSE, height = 12, width = 12)
ggplot(osteoblast, aes(x = enrichment, y = -log10(p.adjust)))+
geom_point(aes(size=observed), color = "grey")+
theme_bw()+
ylim(0,4)+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    text = element_text(size = 20))+
    geom_label_repel(data = subset(osteoblast, rank(osteoblast$p.adjust) <= 15), aes(label = Pathway), size =5)
dev.off()



