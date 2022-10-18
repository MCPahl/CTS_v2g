library(tidyverse)
library(GenomicRanges)

#Load SNP data
leadSNPs = read.delim("Skuladottir2022_carpaltunnel_leadSNPs.txt") %>%
 dplyr::rename(sentinel = rs_name)

proxies = read.delim("combined_query_snp_list.txt", row.names=NULL)
proxies = proxies %>% filter(R2 > 0.8) %>% 
 dplyr::rename(sentinel = query_snp, proxy = RS_Number, proxy_coord_hg19 = Coord, alleles = Alleles, maf = MAF, sentinel_proxy_distance = Distance)

proxies.gr = GRanges(seqnames= gsub(":.*", "", proxies$proxy_coord_hg19), IRanges(as.integer(gsub(".*:", "", proxies$proxy_coord_hg19)), as.integer(gsub(".*:", "", proxies$proxy_coord_hg19))))

db_dir = "/mnt/isilon/sfgi/shinyapp_GCI_db/db"

db = list.files(db_dir)
query_cells = unique(gsub("\\..*", "", db))
query_cells = query_cells[query_cells %in% c("name_convert", "initial_load")==F]


query_cell = "PHMB_Myotube"
v2g_all = list()[1:length(query_cells)]
v2g_all = lapply(query_cells, function(query_cell){
	files = db[grepl(paste0(query_cell, "\\."), db)]
	paths = paste0(db_dir, "/", files)

	atac = read.delim(paths[1])
	frag_int = read.delim(paths[2])
	gene_frag = read.delim(paths[3])
	open_frag = read.delim(paths[4])
	open_promoter = read.delim(paths[5])

	open_interactions_index = as.data.frame(findOverlaps(GRanges(open_frag), proxies.gr))
	proxies_open_frag = data.frame(proxies[open_interactions_index[,2],], open_frag[open_interactions_index[,1],])
	#a = proxies_open_frag %>% dplyr::rename(bait_id = frag_id) %>% left_join(gene_frag) %>% mutate(snp_pos = "bait")
	b = proxies_open_frag %>% dplyr::rename(oe_id = frag_id) %>% left_join(gene_frag) %>% mutate(snp_pos = "oe")

	v2g = rbind(b) %>% select(sentinel, proxy, proxy_coord_hg19,  alleles, maf, sentinel_proxy_distance, chr, start,end, id, frag, cell, gene_name, gene_id, biotype,  oeGenes, snp_pos)
	v2g = v2g[(is.na(v2g$gene_name) & is.na(v2g$oeGenes))==F,]
	v2g = v2g[(v2g$snp_pos=="bait" &  is.na(v2g$oeGenes))==F,]
	v2g
})
v2g_all = do.call("rbind", v2g_all)

v2g_all = left_join(v2g_all, leadSNPs)%>% 
 select(-Chr.Position_hg38) %>% 
 relocate(sentinel,Closest_PC, CodingEffect,  CodingChange, EA, OA, sentinel, proxy, proxy_coord_hg19,alleles, maf, sentinel_proxy_distance, chr, start, end, id, frag, cell, gene_name, gene_id, biotype, oeGenes, snp_pos)



write.table(unique(v2g_all), file = "v2g_carpal_tunnel_sfg_database.txt", quote=F, row.names=F, sep="\t")

pdf("v2g_gene_counts.pdf")
celltype_gene_count = v2g_all %>% select(cell, gene_id) %>% group_by(cell) %>% unique() %>% tally()
ggplot(celltype_gene_count, aes(x= cell, y = n))+
geom_bar(stat="identity")+
coord_flip()+
theme_bw()
dev.off()
