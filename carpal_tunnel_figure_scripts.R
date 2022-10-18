library(tidyverse)
v2g = read.delim("v2g_carpal_tunnel_sfg_database.txt")


celltypes = c("PHMB_Myotube", "hMSC_BMP2")
v2g = v2g %>% 
	filter(cell %in% celltypes)


#Prep how many sentinels contact nearest gene
t = v2g %>% 
	select(sentinel, Closest_PC, cell, gene_name) %>%
	distinct() %>%
	mutate(class = ifelse(Closest_PC == gene_name, "closest", "not_closest")) %>% 
	group_by(sentinel, cell) %>%
	summarize(class = ifelse(Reduce("|", class %in% "closest") & length(unique(.)) == 1, "only_closest", ifelse(Reduce("|", class %in% "closest"), "includes_closest", "excludes_closest"))) %>%
	  as.data.frame()


ggplot()


Stim_HIPK1 = STIM