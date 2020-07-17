matrix_enzymes = read.table(file.choose(), header = TRUE, sep = ";", row.names = 1)
mt_enzymes_scaled = as.matrix(scale(matrix_enzymes))
install.packages("gplots")
library("gplots")
pdf("Documents/proyecto_tesis/analisis_rutas/dendrogramas.pdf")
par(mfrow=c(6,4),mar=c(1,1,1,1), oma=c(0,0,0,12))
heatmap.2(mt_enzymes_scaled, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")

par(mfrow=c(6,4),mar=c(1,1,1,1), oma=c(0,0,0,12))
heatmap.2(mt_enzymes_scaled, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram = "column")

par(mfrow=c(6,4),mar=c(1,1,1,1), oma=c(0,0,0,12))
heatmap.2(mt_enzymes_scaled, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram = "row")
dev.off()