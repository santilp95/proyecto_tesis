datos_genes_sin_filtro<-read.table(file.choose(), header =  TRUE, sep=",")
datos_genes_con_filtro<-read.csv(header = TRUE, file.choose())

#Filtro con valor p <0.05 y log2Foldchange >2#

AvrBs2_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$AvrBs2.Control___padj <= 0.05 )
Up_AvrBs2_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$AvrBs2.Control___padj <= 0.05 & datos_genes_sin_filtro$AvrBs2.Control___log2FoldChange >= 2)
Down_AvrBs2_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$AvrBs2.Control___padj <= 0.05 & datos_genes_sin_filtro$AvrBs2.Control___log2FoldChange <= -2)

XB3.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XB3.Control___padj <= 0.05)
Up_XB3_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XB3.Control___padj <= 0.05 & datos_genes_sin_filtro$XB3.Control___log2FoldChange >= 2)
Down_XB3_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XB3.Control___padj <= 0.05 & datos_genes_sin_filtro$XB3.Control___log2FoldChange <= -2)

XopAE.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAE.Control___padj <= 0.05)
Up_XopAE_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAE.Control___padj <= 0.05 & datos_genes_sin_filtro$XopAE.Control___log2FoldChange >= 2)
Down_XopAE_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAE.Control___padj <= 0.05 & datos_genes_sin_filtro$XopAE.Control___log2FoldChange <= -2)

XopAG.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAG.Control___padj <= 0.05)
Up_XopAG_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAG.Control___padj <= 0.05 & datos_genes_sin_filtro$XopAG.Control___log2FoldChange >= 2)
Down_XopAG_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAG.Control___padj <= 0.05 & datos_genes_sin_filtro$XopAG.Control___log2FoldChange <= -2)

XopAK.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAK.Control___padj <= 0.05)
Up_XopAK_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAK.Control___padj <= 0.05 & datos_genes_sin_filtro$XopAK.Control___log2FoldChange >= 2)
Down_XopAK_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAK.Control___padj <= 0.05 & datos_genes_sin_filtro$XopAK.Control___log2FoldChange <= -2)

XopAO1.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAO1.Control___padj <= 0.05)
Up_XopAO1_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAO1.Control___padj <= 0.05 & datos_genes_sin_filtro$XopAO1.Control___log2FoldChange >= 2)
Down_XopAO1_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAO1.Control___padj <= 0.05 & datos_genes_sin_filtro$XopAO1.Control___log2FoldChange <= -2)

XopC2.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopC2.Control___padj <= 0.05)
Up_XopC2_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopC2.Control___padj <= 0.05 & datos_genes_sin_filtro$XopC2.Control___log2FoldChange >= 2)
Down_XopC2_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopC2.Control___padj <= 0.05 & datos_genes_sin_filtro$XopC2.Control___log2FoldChange <= -2)

XopE1.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE1.Control___padj <= 0.05)
Up_XopE1_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE1.Control___padj <= 0.05 & datos_genes_sin_filtro$XopE1.Control___log2FoldChange >= 2)
Down_XopE1_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE1.Control___padj <= 0.05 & datos_genes_sin_filtro$XopE1.Control___log2FoldChange <= -2)

XopE4.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE4.Control___padj <= 0.05)
Up_XopE4_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE4.Control___padj <= 0.05 & datos_genes_sin_filtro$XopE4.Control___log2FoldChange >= 2)
Down_XopE4_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE4.Control___padj <= 0.05 & datos_genes_sin_filtro$XopE4.Control___log2FoldChange <= -2)

XopK.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopK.Control___padj <= 0.05)
Up_XopK_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopK.Control___padj <= 0.05 & datos_genes_sin_filtro$XopK.Control___log2FoldChange >= 2)
Down_XopK_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopK.Control___padj <= 0.05 & datos_genes_sin_filtro$XopK.Control___log2FoldChange <= -2)

XopL.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopL.Control___padj <= 0.05)
Up_XopL_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopL.Control___padj <= 0.05 & datos_genes_sin_filtro$XopL.Control___log2FoldChange >= 2)
Down_XopL_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopL.Control___padj <= 0.05 & datos_genes_sin_filtro$XopL.Control___log2FoldChange <= -2)

XopN.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopN.Control___padj <= 0.05)
Up_XopN_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopN.Control___padj <= 0.05 & datos_genes_sin_filtro$XopN.Control___log2FoldChange >= 2)
Down_XopN_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopN.Control___padj <= 0.05 & datos_genes_sin_filtro$XopN.Control___log2FoldChange <= -2)

XopQ.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopQ.Control___padj <= 0.05)
Up_XopQ_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopQ.Control___padj <= 0.05 & datos_genes_sin_filtro$XopQ.Control___log2FoldChange >= 2)
Down_XopQ_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopQ.Control___padj <= 0.05 & datos_genes_sin_filtro$XopQ.Control___log2FoldChange <= -2)

XopR.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopR.Control___padj <= 0.05)
Up_XopR_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopR.Control___padj <= 0.05 & datos_genes_sin_filtro$XopR.Control___log2FoldChange >= 2)
Down_XopR_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopR.Control___padj <= 0.05 & datos_genes_sin_filtro$XopR.Control___log2FoldChange <= -2)

XopV.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopV.Control___padj <= 0.05)
Up_XopV_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopV.Control___padj <= 0.05 & datos_genes_sin_filtro$XopV.Control___log2FoldChange >= 2)
Down_XopV_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopV.Control___padj <= 0.05 & datos_genes_sin_filtro$XopV.Control___log2FoldChange <= -2)

XopZ.Control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopZ.Control___padj <= 0.05)
Up_XopZ_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopZ.Control___padj <= 0.05 & datos_genes_sin_filtro$XopZ.Control___log2FoldChange >= 2)
Down_XopZ_control_p_sub<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopZ.Control___padj <= 0.05 & datos_genes_sin_filtro$XopZ.Control___log2FoldChange <= -2)


#Filtro con valor p <0.01 y log2Foldchange >2#

AvrBs2_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$AvrBs2.Control___padj <= 0.01 )
Up_AvrBs2_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$AvrBs2.Control___padj <= 0.01 & datos_genes_sin_filtro$AvrBs2.Control___log2FoldChange >= 2)
Down_AvrBs2_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$AvrBs2.Control___padj <= 0.01 & datos_genes_sin_filtro$AvrBs2.Control___log2FoldChange <= -2)

XB3.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XB3.Control___padj <= 0.01)
Up_XB3_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XB3.Control___padj <= 0.01 & datos_genes_sin_filtro$XB3.Control___log2FoldChange >= 2)
Down_XB3_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XB3.Control___padj <= 0.01 & datos_genes_sin_filtro$XB3.Control___log2FoldChange <= -2)

XopAE.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAE.Control___padj <= 0.01)
Up_XopAE_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAE.Control___padj <= 0.01 & datos_genes_sin_filtro$XopAE.Control___log2FoldChange >= 2)
Down_XopAE_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAE.Control___padj <= 0.01 & datos_genes_sin_filtro$XopAE.Control___log2FoldChange <= -2)

XopAG.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAG.Control___padj <= 0.01)
Up_XopAG_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAG.Control___padj <= 0.01 & datos_genes_sin_filtro$XopAG.Control___log2FoldChange >= 2)
Down_XopAG_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAG.Control___padj <= 0.01 & datos_genes_sin_filtro$XopAG.Control___log2FoldChange <= -2)

XopAK.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAK.Control___padj <= 0.01)
Up_XopAK_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAK.Control___padj <= 0.01 & datos_genes_sin_filtro$XopAK.Control___log2FoldChange >= 2)
Down_XopAK_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAK.Control___padj <= 0.01 & datos_genes_sin_filtro$XopAK.Control___log2FoldChange <= -2)

XopAO1.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAO1.Control___padj <= 0.01)
Up_XopAO1_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAO1.Control___padj <= 0.01 & datos_genes_sin_filtro$XopAO1.Control___log2FoldChange >= 2)
Down_XopAO1_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopAO1.Control___padj <= 0.01 & datos_genes_sin_filtro$XopAO1.Control___log2FoldChange <= -2)

XopC2.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopC2.Control___padj <= 0.01)
Up_XopC2_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopC2.Control___padj <= 0.01 & datos_genes_sin_filtro$XopC2.Control___log2FoldChange >= 2)
Down_XopC2_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopC2.Control___padj <= 0.01 & datos_genes_sin_filtro$XopC2.Control___log2FoldChange <= -2)

XopE1.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE1.Control___padj <= 0.01)
Up_XopE1_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE1.Control___padj <= 0.01 & datos_genes_sin_filtro$XopE1.Control___log2FoldChange >= 2)
Down_XopE1_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE1.Control___padj <= 0.01 & datos_genes_sin_filtro$XopE1.Control___log2FoldChange <= -2)

XopE4.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE4.Control___padj <= 0.01)
Up_XopE4_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE4.Control___padj <= 0.01 & datos_genes_sin_filtro$XopE4.Control___log2FoldChange >= 2)
Down_XopE4_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopE4.Control___padj <= 0.01 & datos_genes_sin_filtro$XopE4.Control___log2FoldChange <= -2)

XopK.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopK.Control___padj <= 0.01)
Up_XopK_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopK.Control___padj <= 0.01 & datos_genes_sin_filtro$XopK.Control___log2FoldChange >= 2)
Down_XopK_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopK.Control___padj <= 0.01 & datos_genes_sin_filtro$XopK.Control___log2FoldChange <= -2)

XopL.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopL.Control___padj <= 0.01)
Up_XopL_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopL.Control___padj <= 0.01 & datos_genes_sin_filtro$XopL.Control___log2FoldChange >= 2)
Down_XopL_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopL.Control___padj <= 0.01 & datos_genes_sin_filtro$XopL.Control___log2FoldChange <= -2)

XopN.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopN.Control___padj <= 0.01)
Up_XopN_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopN.Control___padj <= 0.01 & datos_genes_sin_filtro$XopN.Control___log2FoldChange >= 2)
Down_XopN_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopN.Control___padj <= 0.01 & datos_genes_sin_filtro$XopN.Control___log2FoldChange <= -2)

XopQ.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopQ.Control___padj <= 0.01)
Up_XopQ_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopQ.Control___padj <= 0.01 & datos_genes_sin_filtro$XopQ.Control___log2FoldChange >= 2)
Down_XopQ_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopQ.Control___padj <= 0.01 & datos_genes_sin_filtro$XopQ.Control___log2FoldChange <= -2)

XopR.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopR.Control___padj <= 0.01)
Up_XopR_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopR.Control___padj <= 0.01 & datos_genes_sin_filtro$XopR.Control___log2FoldChange >= 2)
Down_XopR_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopR.Control___padj <= 0.01 & datos_genes_sin_filtro$XopR.Control___log2FoldChange <= -2)

XopV.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopV.Control___padj <= 0.01)
Up_XopV_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopV.Control___padj <= 0.01 & datos_genes_sin_filtro$XopV.Control___log2FoldChange >= 2)
Down_XopV_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopV.Control___padj <= 0.01 & datos_genes_sin_filtro$XopV.Control___log2FoldChange <= -2)

XopZ.Control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopZ.Control___padj <= 0.01)
Up_XopZ_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopZ.Control___padj <= 0.01 & datos_genes_sin_filtro$XopZ.Control___log2FoldChange >= 2)
Down_XopZ_control_p_sub_0.01<-subset(datos_genes_sin_filtro$Row.names, datos_genes_sin_filtro$XopZ.Control___padj <= 0.01 & datos_genes_sin_filtro$XopZ.Control___log2FoldChange <= -2)

#Listados de genes en orden descendiente por LOG_FOLD_CHANGE >2#

library(ggplot2)
library(dplyr)

#AvrBs2.Control#
#Esto se utiliza para filtrar todo el data frame, como un subset pero mantiene todas las variables#
AvrBs2.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, AvrBs2.Control___padj <= 0.01 & AvrBs2.Control___log2FoldChange >= 2)
#keeps es un vector que identifica variables para dejar, y luego con los corchetes se utiliza este vector, posterior al df, para dejarlas#
keeps1=c("Row.names", "AvrBs2.Control___padj", "AvrBs2.Control___log2FoldChange")
AvrBs2.Control___padj_0.01_LOG_2_filt=AvrBs2.Control___padj_0.01_LOG_2[keeps1]
AvrBs2.Control___padj_0.01_LOG_2_filt
AvrBs2.Control___padj_0.01_LOG_2_filt_desc=arrange(AvrBs2.Control___padj_0.01_LOG_2_filt, desc(AvrBs2.Control___log2FoldChange))
write.table(AvrBs2.Control___padj_0.01_LOG_2_filt_desc, file = "AvrBs2 list", sep = ",", row.names = FALSE)
quantile(AvrBs2.Control___padj_0.01_LOG_2_filt_desc$AvrBs2.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data = AvrBs2.Control___padj_0.01_LOG_2_filt_desc, aes(x = AvrBs2.Control___log2FoldChange)) +
  geom_histogram()

datos_genes_sin_filtro

#XB3.Control#
XB3.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XB3.Control___padj <= 0.01 & XB3.Control___log2FoldChange >= 2)
keeps2=c("Row.names", "XB3.Control___padj", "XB3.Control___log2FoldChange")
XB3.Control___padj_0.01_LOG_2_filt=XB3.Control___padj_0.01_LOG_2[keeps2]
XB3.Control___padj_0.01_LOG_2_filt
XB3.Control___padj_0.01_LOG_2_filt_desc=arrange(XB3.Control___padj_0.01_LOG_2_filt, desc(XB3.Control___log2FoldChange))
write.table(XB3.Control___padj_0.01_LOG_2_filt_desc, file = "XB3 list", sep = ",", row.names = FALSE)
quantile(XB3.Control___padj_0.01_LOG_2_filt_desc$XB3.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data = XB3.Control___padj_0.01_LOG_2_filt_desc, aes(x = XB3.Control___log2FoldChange)) +
  geom_histogram()


"XopAE.Control"

XopAE.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopAE.Control___padj <= 0.01 & XopAE.Control___log2FoldChange >= 2)
keeps3=c("Row.names", "XopAE.Control___padj", "XopAE.Control___log2FoldChange")
XopAE.Control___padj_0.01_LOG_2_filt=XopAE.Control___padj_0.01_LOG_2[keeps3]
XopAE.Control___padj_0.01_LOG_2_filt
XopAE.Control___padj_0.01_LOG_2_filt_desc=arrange(XopAE.Control___padj_0.01_LOG_2_filt, desc(XopAE.Control___log2FoldChange))
write.table(XopAE.Control___padj_0.01_LOG_2_filt_desc, file = "XopAE list", sep = ",", row.names = FALSE)
quantile(XopAE.Control___padj_0.01_LOG_2_filt_desc$XopAE.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data = XopAE.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopAE.Control___log2FoldChange)) +
  geom_histogram()

#XopAG.Control#

XopAG.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopAG.Control___padj <= 0.01 & XopAG.Control___log2FoldChange >= 2)
keeps4=c("Row.names", "XopAG.Control___padj", "XopAG.Control___log2FoldChange")
XopAG.Control___padj_0.01_LOG_2_filt=XopAG.Control___padj_0.01_LOG_2[keeps4]
XopAG.Control___padj_0.01_LOG_2_filt
XopAG.Control___padj_0.01_LOG_2_filt_desc=arrange(XopAG.Control___padj_0.01_LOG_2_filt, desc(XopAG.Control___log2FoldChange))
write.table(XopAG.Control___padj_0.01_LOG_2_filt_desc, file = "XopAG list", sep = ",", row.names = FALSE)
quantile(XopAG.Control___padj_0.01_LOG_2_filt_desc$XopAG.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data = XopAG.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopAG.Control___log2FoldChange)) +
  geom_histogram()

#XopAK.Control#

XopAK.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopAK.Control___padj <= 0.01 & XopAK.Control___log2FoldChange >= 2)
keeps5=c("Row.names", "XopAK.Control___padj", "XopAK.Control___log2FoldChange")
XopAK.Control___padj_0.01_LOG_2_filt=XopAK.Control___padj_0.01_LOG_2[keeps5]
XopAK.Control___padj_0.01_LOG_2_filt
XopAK.Control___padj_0.01_LOG_2_filt_desc=arrange(XopAK.Control___padj_0.01_LOG_2_filt, desc(XopAK.Control___log2FoldChange))
write.table(XopAK.Control___padj_0.01_LOG_2_filt_desc, file = "XopAK list", sep = ",", row.names = FALSE)
quantile(XopAK.Control___padj_0.01_LOG_2_filt_desc$XopAK.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data = XopAK.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopAK.Control___log2FoldChange)) +
  geom_histogram()

#XopAO1#

XopAO1.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopAO1.Control___padj <= 0.01 & XopAO1.Control___log2FoldChange >= 2)
keeps6=c("Row.names", "XopAO1.Control___padj", "XopAO1.Control___log2FoldChange")
XopAO1.Control___padj_0.01_LOG_2_filt=XopAO1.Control___padj_0.01_LOG_2[keeps6]
XopAO1.Control___padj_0.01_LOG_2_filt
XopAO1.Control___padj_0.01_LOG_2_filt_desc=arrange(XopAO1.Control___padj_0.01_LOG_2_filt, desc(XopAO1.Control___log2FoldChange))
write.table(XopAO1.Control___padj_0.01_LOG_2_filt_desc, file = "XopAO1 list", sep = ",", row.names = FALSE)
quantile(XopAO1.Control___padj_0.01_LOG_2_filt_desc$XopAO1.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopAO1.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopAO1.Control___log2FoldChange)) +
  geom_histogram()

#XopC2#

XopC2.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopC2.Control___padj <= 0.01 & XopC2.Control___log2FoldChange >= 2)
keeps7=c("Row.names", "XopC2.Control___padj", "XopC2.Control___log2FoldChange")
XopC2.Control___padj_0.01_LOG_2_filt=XopC2.Control___padj_0.01_LOG_2[keeps7]
XopC2.Control___padj_0.01_LOG_2_filt
XopC2.Control___padj_0.01_LOG_2_filt_desc=arrange(XopC2.Control___padj_0.01_LOG_2_filt, desc(XopC2.Control___log2FoldChange))
write.table(XopC2.Control___padj_0.01_LOG_2_filt_desc, file = "XopC2 list", sep = ",", row.names = FALSE)
quantile(XopC2.Control___padj_0.01_LOG_2_filt_desc$XopC2.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopC2.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopC2.Control___log2FoldChange)) +
  geom_histogram()

#XopE1#

XopE1.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopE1.Control___padj <= 0.01 & XopE1.Control___log2FoldChange >= 2)
keeps8=c("Row.names", "XopE1.Control___padj", "XopE1.Control___log2FoldChange")
XopE1.Control___padj_0.01_LOG_2_filt=XopE1.Control___padj_0.01_LOG_2[keeps8]
XopE1.Control___padj_0.01_LOG_2_filt
XopE1.Control___padj_0.01_LOG_2_filt_desc=arrange(XopE1.Control___padj_0.01_LOG_2_filt, desc(XopE1.Control___log2FoldChange))
write.table(XopE1.Control___padj_0.01_LOG_2_filt_desc, file = "XopE1 list", sep = ",", row.names = FALSE)
quantile(XopE1.Control___padj_0.01_LOG_2_filt_desc$XopE1.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopE1.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopE1.Control___log2FoldChange)) +
  geom_histogram()

#XopE4#

XopE4.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopE4.Control___padj <= 0.01 & XopE4.Control___log2FoldChange >= 2)
keeps9=c("Row.names", "XopE4.Control___padj", "XopE4.Control___log2FoldChange")
XopE4.Control___padj_0.01_LOG_2_filt=XopE4.Control___padj_0.01_LOG_2[keeps9]
XopE4.Control___padj_0.01_LOG_2_filt
XopE4.Control___padj_0.01_LOG_2_filt_desc=arrange(XopE4.Control___padj_0.01_LOG_2_filt, desc(XopE4.Control___log2FoldChange))
write.table(XopE4.Control___padj_0.01_LOG_2_filt_desc, file = "XopE4 list", sep = ",", row.names = FALSE)
quantile(XopE4.Control___padj_0.01_LOG_2_filt_desc$XopE4.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopE4.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopE4.Control___log2FoldChange)) +
  geom_histogram()


#XopK#

XopK.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopK.Control___padj <= 0.01 & XopK.Control___log2FoldChange >= 2)
keeps10=c("Row.names", "XopK.Control___padj", "XopK.Control___log2FoldChange")
XopK.Control___padj_0.01_LOG_2_filt=XopK.Control___padj_0.01_LOG_2[keeps10]
XopK.Control___padj_0.01_LOG_2_filt
XopK.Control___padj_0.01_LOG_2_filt_desc=arrange(XopK.Control___padj_0.01_LOG_2_filt, desc(XopK.Control___log2FoldChange))
write.table(XopK.Control___padj_0.01_LOG_2_filt_desc, file = "XopK list", sep = ",", row.names = FALSE)
quantile(XopK.Control___padj_0.01_LOG_2_filt_desc$XopK.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopK.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopK.Control___log2FoldChange)) +
  geom_histogram()


#XopL#

XopL.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopL.Control___padj <= 0.01 & XopL.Control___log2FoldChange >= 2)
keeps11=c("Row.names", "XopL.Control___padj", "XopL.Control___log2FoldChange")
XopL.Control___padj_0.01_LOG_2_filt=XopL.Control___padj_0.01_LOG_2[keeps11]
XopL.Control___padj_0.01_LOG_2_filt
XopL.Control___padj_0.01_LOG_2_filt_desc=arrange(XopL.Control___padj_0.01_LOG_2_filt, desc(XopL.Control___log2FoldChange))
write.table(XopL.Control___padj_0.01_LOG_2_filt_desc, file = "XopL list", sep = ",", row.names = FALSE)
quantile(XopL.Control___padj_0.01_LOG_2_filt_desc$XopL.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopL.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopL.Control___log2FoldChange)) +
  geom_histogram()

#XopN#

XopN.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopN.Control___padj <= 0.01 & XopN.Control___log2FoldChange >= 2)
keeps12=c("Row.names", "XopN.Control___padj", "XopN.Control___log2FoldChange")
XopN.Control___padj_0.01_LOG_2_filt=XopN.Control___padj_0.01_LOG_2[keeps12]
XopN.Control___padj_0.01_LOG_2_filt
XopN.Control___padj_0.01_LOG_2_filt_desc=arrange(XopN.Control___padj_0.01_LOG_2_filt, desc(XopN.Control___log2FoldChange))
write.table(XopN.Control___padj_0.01_LOG_2_filt_desc, file = "XopN list", sep = ",", row.names = FALSE)
quantile(XopN.Control___padj_0.01_LOG_2_filt_desc$XopN.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopN.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopN.Control___log2FoldChange)) +
  geom_histogram()

#XopQ#

XopQ.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopQ.Control___padj <= 0.01 & XopQ.Control___log2FoldChange >= 2)
keeps13=c("Row.names", "XopQ.Control___padj", "XopQ.Control___log2FoldChange")
XopQ.Control___padj_0.01_LOG_2_filt=XopQ.Control___padj_0.01_LOG_2[keeps13]
XopQ.Control___padj_0.01_LOG_2_filt
XopQ.Control___padj_0.01_LOG_2_filt_desc=arrange(XopQ.Control___padj_0.01_LOG_2_filt, desc(XopQ.Control___log2FoldChange))
write.table(XopQ.Control___padj_0.01_LOG_2_filt_desc, file = "XopQ list", sep = ",", row.names = FALSE)
quantile(XopQ.Control___padj_0.01_LOG_2_filt_desc$XopQ.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopQ.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopQ.Control___log2FoldChange)) +
  geom_histogram()

#XopR#

XopR.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopR.Control___padj <= 0.01 & XopR.Control___log2FoldChange >= 2)
keeps14=c("Row.names", "XopR.Control___padj", "XopR.Control___log2FoldChange")
XopR.Control___padj_0.01_LOG_2_filt=XopR.Control___padj_0.01_LOG_2[keeps14]
XopR.Control___padj_0.01_LOG_2_filt
XopR.Control___padj_0.01_LOG_2_filt_desc=arrange(XopR.Control___padj_0.01_LOG_2_filt, desc(XopR.Control___log2FoldChange))
write.table(XopR.Control___padj_0.01_LOG_2_filt_desc, file = "XopR list", sep = ",", row.names = FALSE)
quantile(XopR.Control___padj_0.01_LOG_2_filt_desc$XopR.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopR.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopR.Control___log2FoldChange)) +
  geom_histogram()

#XopV#

XopV.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopV.Control___padj <= 0.01 & XopV.Control___log2FoldChange >= 2)
keeps15=c("Row.names", "XopV.Control___padj", "XopV.Control___log2FoldChange")
XopV.Control___padj_0.01_LOG_2_filt=XopV.Control___padj_0.01_LOG_2[keeps15]
XopV.Control___padj_0.01_LOG_2_filt
XopV.Control___padj_0.01_LOG_2_filt_desc=arrange(XopV.Control___padj_0.01_LOG_2_filt, desc(XopV.Control___log2FoldChange))
write.table(XopV.Control___padj_0.01_LOG_2_filt_desc, file = "XopV list", sep = ",", row.names = FALSE)
quantile(XopV.Control___padj_0.01_LOG_2_filt_desc$XopV.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopV.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopV.Control___log2FoldChange)) +
  geom_histogram()

#XopZ#

XopZ.Control___padj_0.01_LOG_2=filter(datos_genes_sin_filtro, XopZ.Control___padj <= 0.01 & XopZ.Control___log2FoldChange >= 2)
keeps16=c("Row.names", "XopZ.Control___padj", "XopZ.Control___log2FoldChange")
XopZ.Control___padj_0.01_LOG_2_filt=XopZ.Control___padj_0.01_LOG_2[keeps16]
XopZ.Control___padj_0.01_LOG_2_filt
XopZ.Control___padj_0.01_LOG_2_filt_desc=arrange(XopZ.Control___padj_0.01_LOG_2_filt, desc(XopZ.Control___log2FoldChange))
write.table(XopZ.Control___padj_0.01_LOG_2_filt_desc, file = "XopZ list", sep = ",", row.names = FALSE)
quantile(XopZ.Control___padj_0.01_LOG_2_filt_desc$XopZ.Control___log2FoldChange, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

ggplot(data =XopZ.Control___padj_0.01_LOG_2_filt_desc, aes(x = XopZ.Control___log2FoldChange)) +
  geom_histogram()


#Heatmaps#

para_graficar_logfc=datos_genes_sin_filtro
para_graficar_logfc
para_graficar_logfc$Symbol= NULL
para_graficar_logfc$AvrBs2.Control___padj = NULL
para_graficar_logfc$XB3.Control___padj= NULL
para_graficar_logfc$XopAE.Control___padj=NULL
para_graficar_logfc$XopAG.Control___padj=NULL
para_graficar_logfc$XopAK.Control___padj=NULL
para_graficar_logfc$XopAO1.Control___padj=NULL
para_graficar_logfc$XopC2.Control___padj=NULL
para_graficar_logfc$XopE1.Control___padj=NULL
para_graficar_logfc$XopE4.Control___padj=NULL
para_graficar_logfc$XopK.Control___padj=NULL
para_graficar_logfc$XopL.Control___padj=NULL
para_graficar_logfc$XopN.Control___padj=NULL
para_graficar_logfc$XopQ.Control___padj=NULL
para_graficar_logfc$XopR.Control___padj=NULL
para_graficar_logfc$XopV.Control___padj=NULL
para_graficar_logfc$XopZ.Control___padj=NULL
para_graficar_logfc$data=NULL
para_graficar_logfc$Control_1=NULL
para_graficar_logfc$Control_2=NULL
para_graficar_logfc$Control_3=NULL
para_graficar_logfc$XopN_1=NULL
para_graficar_logfc$XopN_2=NULL
para_graficar_logfc$XopN_3=NULL
para_graficar_logfc$XopE4_1=NULL
para_graficar_logfc$XopE4_2=NULL
para_graficar_logfc$XopE4_3=NULL
para_graficar_logfc$XopC2_1=NULL
para_graficar_logfc$XopC2_2=NULL
para_graficar_logfc$XopC2_3=NULL
para_graficar_logfc$XB3_1=NULL
para_graficar_logfc$XB3_2=NULL
para_graficar_logfc$XB3_3=NULL
para_graficar_logfc$XopQ_1=NULL
para_graficar_logfc$XopQ_2=NULL
para_graficar_logfc$XopQ_3=NULL
para_graficar_logfc$XopL_1=NULL
para_graficar_logfc$XopL_2=NULL
para_graficar_logfc$XopL_3=NULL
para_graficar_logfc$XopAG_1=NULL
para_graficar_logfc$XopAG_2=NULL
para_graficar_logfc$XopAG_3=NULL
para_graficar_logfc$XopE1_1=NULL
para_graficar_logfc$XopE1_2=NULL
para_graficar_logfc$XopE1_3=NULL
para_graficar_logfc$XopK_1=NULL
para_graficar_logfc$XopK_2=NULL
para_graficar_logfc$XopK_3=NULL
para_graficar_logfc$AvrBs2_1=NULL
para_graficar_logfc$AvrBs2_2=NULL
para_graficar_logfc$AvrBs2_3=NULL
para_graficar_logfc$XopV_1=NULL
para_graficar_logfc$XopV_2=NULL
para_graficar_logfc$XopV_3=NULL
para_graficar_logfc$XopR_1=NULL
para_graficar_logfc$XopR_2=NULL
para_graficar_logfc$XopR_3=NULL
para_graficar_logfc$XopAK_1=NULL
para_graficar_logfc$XopAK_2=NULL
para_graficar_logfc$XopAK_3=NULL
para_graficar_logfc$XopAE_1=NULL
para_graficar_logfc$XopAE_2=NULL
para_graficar_logfc$XopAE_3=NULL
para_graficar_logfc$XopAO1_1=NULL
para_graficar_logfc$XopAO2_1=NULL
para_graficar_logfc$XopAO3_1=NULL
para_graficar_logfc$XopAO1_2=NULL
para_graficar_logfc$XopAO1_3=NULL
para_graficar_logfc$XopZ_1=NULL
para_graficar_logfc$XopZ_2=NULL
para_graficar_logfc

y=data.matrix(para_graficar_logfc_d)

para_graficar_logfc_d= para_graficar_logfc[, -1]
rownames(para_graficar_logfc_d) = para_graficar_logfc[, 1]
para_graficar_logfc_d 
is.na(para_graficar_logfc_d) <- sapply(para_graficar_logfc_d, is.infinite)
para_graficar_logfc_d[is.na(para_graficar_logfc_d)] <- 0
para_graficar_logfc_d[is.nan(para_graficar_logfc_d)] <- 0

para_graficar_logfc_d

library("gplots") 


expresion_diferecial_Log2FC=heatmap(y, scale = "none", main = "Discriminaci�n seg�n Log2FoldChange", trace =  "none", margins = c(10,12), cexRow=2)

pheatmap(y, cutree_rows = 4)
