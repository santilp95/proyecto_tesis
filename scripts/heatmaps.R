setwd("Documents/proyecto_tesis/analisis_rutas")
getwd()
matrix_enzymes = read.table("enzimas_rutas_completas.csv", header = TRUE, sep = ",", row.names = 1)
matrix_groups= read.table("cluster_rutas.csv", header = TRUE, sep = ";", row.names = 1)
mt_enzymes_scaled = as.matrix(scale(matrix_enzymes, center = FALSE, scale = apply(matrix_enzymes, 2, sd, na.rm = TRUE)))
mt_groups_scaled = as.matrix(scale(matrix_groups, center = FALSE, scale = apply(matrix_enzymes, 2, sd, na.rm = TRUE)))
library("gplots")
pdf("dendrogramas_rutas_completas.pdf", height = 50)
par(mfrow=c(6,4),mar=c(1,1,1,1), oma=c(0,0,0,12))
heatmap.2(mt_enzymes_scaled, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram = "both")

par(mfrow=c(6,4),mar=c(1,1,1,1), oma=c(0,0,0,12))
heatmap.2(mt_enzymes_scaled, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram = "column", Rowv = FALSE)

par(mfrow=c(6,4),mar=c(1,1,1,1), oma=c(0,0,0,12))
heatmap.2(mt_enzymes_scaled, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram = "row", Colv=FALSE)
dev.off()

pdf("dendrogramas_grupos_metabolicos.pdf")
par(mfrow=c(6,4),mar=c(1,1,1,1), oma=c(0,0,0,12))
heatmap.2(mt_groups_scaled, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram = "both")

par(mfrow=c(6,4),mar=c(1,1,1,1), oma=c(0,0,0,12))
heatmap.2(mt_groups_scaled, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram = "column", Rowv = FALSE)

par(mfrow=c(6,4),mar=c(1,1,1,1), oma=c(0,0,0,12))
heatmap.2(mt_groups_scaled, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram = "row", Colv=FALSE)
dev.off()
