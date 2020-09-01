#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 18:15:02 2020

@author: ngaitan55
"""
import os
import pandas as pd
os.getcwd()
os.chdir("/home/ngaitan55/Documents/proyecto_tesis/redes_metabolicas/archivos_tab")

data_labels = ["df_ko_AvrBs2.csv", "df_ko_XB3.csv", "df_ko_XopAE.csv", "df_ko_XopAG.csv", "df_ko_XopAK.csv",
               "df_ko_XopAO1.csv", "df_ko_XopC2.csv", "df_ko_XopE1.csv", "df_ko_XopE4.csv", "df_ko_XopK.csv",
               "df_ko_XopL.csv", "df_ko_XopN.csv", "df_ko_XopQ.csv", "df_ko_XopR.csv", "df_ko_XopV.csv"]

for i in data_labels:
    df = pd.read_csv(i, sep = "\t")
    ec = df["EC"] != '0'
    df = df[ec]
    df = df[["EC", "log2_fold_change", "valor_p_ajustado"]]
    df.to_csv(i, index = False, sep = "\t")


