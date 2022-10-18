library(pheatmap)

load("/mnt/isilon/sfgi/pahlm/annotationFiles/gtex/v8/gtex_median_tissue_expression_tpm.RData")


x = read.delim("v2g_carpal_tunnel_sfg_database.txt")
x = x[grepl("Myotube",x$cell),]
muscle_genes = unique(x$gene_name)

gtex.tissue.med.mg = gtex.tissue.med[gtex.tissue.med$gene_name %in% muscle_genes,]
row.id = gtex.tissue.med.mg$gene_name
col.id = names(gtex.tissue.med.mg)
gtex.tissue.med.mg$gene_id = NULL
gtex.tissue.med.mg$gene_name = NULL

gtex.tissue.med.mg = as.matrix(gtex.tissue.med.mg)
row.names(gtex.tissue.med.mg)=row.id


gtex.tissue.med.mg = apply(gtex.tissue.med.mg,1, scale)

index = apply(gtex.tissue.med.mg, 2, is.nan)==F
index = apply(index, 2, function(x){
	Reduce("|", x)
})

gtex.tissue.med.mg = gtex.tissue.med.mg[,index]
row.names(gtex.tissue.med.mg) = col.id[-c(1:2)]
#dir.create("plots")

pdf("plots/gtex_median_expression_heatmap.pdf", width=20, height = 20)
pheatmap(t(gtex.tissue.med.mg), show_colnames = T)
dev.off()