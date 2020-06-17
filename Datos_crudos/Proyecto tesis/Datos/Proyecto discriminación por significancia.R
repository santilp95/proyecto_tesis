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

write.table()



