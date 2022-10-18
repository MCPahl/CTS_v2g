library(VennDiagram)

v2g = read.delim("v2g_carpal_tunnel_sfg_database.txt")

osteoblast = v2g[v2g$cell == "hMSC_BMP2",]
myotube = v2g[v2g$cell == "PHMB_Myotube",]

gene_list = list(Osteoblast = unique(osteoblast$gene_name), Myotube = unique(myotube$gene_name))


venn.diagram(gene_list, 
 	filename="plots/venn.diagram.png",
	output=TRUE,

	fill = c("yellow", "red"), 
	lty = 'blank',
  	imagetype="png" ,
    height = 1028 , 
    width = 1028 , 
    resolution = 300,
    compression = "lzw",

    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.fontfamily = "sans")



