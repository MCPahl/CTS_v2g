library(tidyverse)
load("/mnt/isilon/sfgi/shinyapp_GCI_db/PHMB_Myotube/PHMB_Myotube.Rdata")
gene_of_interest = "BZW2"

gene_of_interest_loops = gene_frag %>%
	filter(gene_name %in% gene_of_interest) %>%
	left_join(frag_int) %>% 
	select("bait_chr", "bait_start", "bait_end", "oe_chr", "oe_start", "oe_end", "score") %>%
	arrange("bait_chr", "bait_start", "bait_end", "oe_chr", "oe_start", "oe_end")


gene_of_interest_loops

write.table(gene_of_interest_loops, file = paste0(gene_of_interest, "_promoter_interacting_region.arcs"), quote=F, row.names=F, sep="\t", col.names=F)




[arcs]
file = BZW2_promoter_interacting_region.arcs
height = 1
color = pink
orientation = inverted
links type = arcs
links line width = 1
show data range = no

#Plot
conda activate py37
pyGenomeTracks  --tracks tracks.ini -o BZW2_myotube.pdf --dpi 300 --region chr7:16124000-16964000