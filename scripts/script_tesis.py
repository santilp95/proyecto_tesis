# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 23:32:35 2020

@author: Usuario
"""


import os

import numpy as np 

import pandas as pd

import matplotlib as mpl

import seaborn as sns


import matplotlib.pyplot as plt

import pylab as pl
from pandas import DataFrame
from matplotlib.backends.backend_pdf import PdfPages
os.getcwd()

import sklearn
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

#Obtener el directorio del proyecto#
os.chdir("/home/ngaitan55/Documents/proyecto_tesis")

#Obtener el ec de clasificación primordial de tipo enzimático a partir de listas de ec#
def obtener_primer_ec(enzimas):
    lista_primer_ec = []
    for x in enzimas:
        lista_primer_ec.append(x[0])
        coma = x.find(',')
        if coma != -1:
            while coma < len(x):
                lista_primer_ec.append(x[coma +1])
                coma = x.find(',', coma +1)
                if coma == -1:
                    break
    return lista_primer_ec
        

#A partir del primer ec se obtiene el tipo de enzima como un str#
def asignar_tipo_enzima(primeros_ec):
    tipos_de_enzimas = ["oxidorreductasa", "transferasa", "hidrolasa", "liasa", "isomerasa", "ligasa", "translocasa"]
    tipos_enzimas_efector = []
    for x in primeros_ec:
        if x == '1':
            tipos_enzimas_efector.append(tipos_de_enzimas[0])
        if x == '2':
            tipos_enzimas_efector.append(tipos_de_enzimas[1])
        if x == '3':
            tipos_enzimas_efector.append(tipos_de_enzimas[2])
        if x == '4':
            tipos_enzimas_efector.append(tipos_de_enzimas[3])
        if x == '5':
            tipos_enzimas_efector.append(tipos_de_enzimas[4])
        if x == '6':
            tipos_enzimas_efector.append(tipos_de_enzimas[5])
        if x == '7':
            tipos_enzimas_efector.append(tipos_de_enzimas[6])
    return tipos_enzimas_efector

#conteo de los tipos enzimáticos según el primer ec#
def contar_tipos_enzimas(lista_primer_ec):
    oxidorreductasas = 0
    transferasas = 0
    hidrolasas = 0
    liasas = 0
    isomerasas = 0
    ligasas = 0
    translocasas = 0
    for x in lista_primer_ec:
        if x == '1':
            oxidorreductasas += 1
        if x == '2':
            transferasas += 1
        if x == '3':
            hidrolasas += 1
        if x == '4':
            liasas += 1
        if x == '5':
            isomerasas += 1
        if x == '6':
            ligasas += 1
        if x == '7':
            translocasas += 1
    df_conteo_tipos_enzimas = [oxidorreductasas, transferasas, hidrolasas, 
                               liasas, isomerasas, ligasas, translocasas]
    return df_conteo_tipos_enzimas

#subir cada query como un df, cambiar NaN por 0 y guardar como un nuevo texto
def cambiar_query(ubicacion_query_interactor: str):
    nuevo_query = pd.read_csv(ubicacion_query_interactor, sep = '\t', header = None)
    nuevo_query = nuevo_query.fillna(0)
    nuevo_query.to_csv(ubicacion_query_interactor, sep = '\t', index = False, header = False)
    return

    
def crear_data_frame_manid_ko_araid_ec(interactor_ID_RAW, query_interactor, no_encontrados):
    KO = []
    if no_encontrados.size != 0:
        if no_encontrados.size != 1:
            for x in no_encontrados:
                index = interactor_ID_RAW[interactor_ID_RAW[0] == x].index
                interactor_ID_RAW.drop(index, inplace = True)
        else:
            x = no_encontrados
            index = interactor_ID_RAW[interactor_ID_RAW[0] == x].index
            interactor_ID_RAW.drop(index, inplace = True)
    for x in interactor_ID_RAW[0]:
        ko = ''
        for y, z in zip(query_interactor[0], query_interactor[1]):
            if x == y:
                if z != '0':
                    ko = z
                else:
                    ko = '0'
                break
        KO.append(ko)
    if len(KO) == len(interactor_ID_RAW):
        df_manihot = pd.DataFrame({'id_manihot': interactor_ID_RAW[0], 'id_arabidopsis': interactor_ID_RAW[1],
                               'EC': interactor_ID_RAW[2], 'KO': KO})
    else:
        return "Hay algún problema con el código"
    return df_manihot

def crossref(lista_buscada: pd.DataFrame, base_datos: pd.DataFrame, ncol_buscada: str, ncol_extraccion: str, unir: bool) -> pd.DataFrame:
    """
    Método que permite hacer un crossref entre una lista con determinados elementos y una base de datos.

    Parameters
    ----------
    lista_buscada : list
        Lista con los elementos que se pretenden buscar en la base de datos.
    base_datos : pd.DataFrame
        Base de datos donde se buscarán los elementos de la lista_buscada.
    ncol_buscada : int
        Dentro de la base de datos, es el indicador de la columna en la cual se buscarán los elmentos de la lista buscada.
    ncol_extraccion : int
        Dentro de la base de datos, es el indicador de la columna de la cual se extraerán los elementos que correspondan
        al atributo por extraer del data point donde el elemento x de la lista buscada es igual (estrictamente) al elemento_bd
        del atributo buscado en la base de datos.
    unir : bool
        True: Si se quiere obtener la base de datos como un inner join
        False: SI se quieren obtener los elementos buscados con su correspondiente elemento encontrado

    Returns
    unir == True: retorna un data frame con sólo los data points que coincidan con la búsqueda. 
    unir == False: retorna un dataframe con la lista buscada y la lista extraida de la base de datos.
    -------

    """
    base_encontrada = pd.merge(lista_buscada, base_datos, how = 'inner', on = ncol_buscada)
    if unir:
        return base_encontrada
    else:
        return base_encontrada[[ncol_buscada, ncol_extraccion]]
    
    
def leer_ec(ubicacion_archivo):
    lista_ecs = []
    with open(ubicacion_archivo, 'r') as lector:
        for linea in lector:
            lista_ecs.append(linea)
    return lista_ecs

def añadir_prefijo_comas(objetivo: list, nuevo_prefijo: str)-> list:
    prefijo= []
    for x in objetivo:
        prefijo.append(nuevo_prefijo)
    objetivo_completo = []
    for x, y in zip(prefijo, objetivo):
        coma_index = y.find(',')
        if coma_index != -1:
            prefijo = x + y[0:coma_index+1]
            resto = ""
            while coma_index < len(y):
                nueva_coma = y.find(',', coma_index+1)
                if nueva_coma != -1:
                    resto += x + y[coma_index+1:nueva_coma+1]
                    coma_index = nueva_coma
                else: 
                    resto += x + y[coma_index+1:len(y)]
                    break
            y = prefijo + resto
            objetivo_completo.append(y)
        else:
            objetivo_completo.append(x+y)
    return objetivo_completo
   
def comparar_valores(lista1: list, lista2: list)-> pd.DataFrame:
    """
    Compara los valores de dos listas y retorna una lista con valores booleanos que indican si son iguales o no

    Parameters
    ----------
    lista1 : list
    lista2 : list

    Returns
    -------
    df : Contiene la lista1, lista2 y su respectiva comparación.

    """
    comparacion = []
    for x, y in zip(lista1, lista2):
        if x == y:
            comparacion.append(True)
        else:
            comparacion.append(False)
    df = pd.concat([lista1, lista2, comparacion])
    return df
            
#Importe de secuencias EC#
AvrBs2 = np.genfromtxt("IDs/AvrBs2/AvrBs2_KEGG.txt", dtype = 'str')
XB3 = np.genfromtxt("IDs/XB3/XB3_KEGG.txt", dtype = 'str')
XopAE = np.genfromtxt("IDs/XopAE/XopAE_KEGG.txt", dtype = 'str')
XopAG = np.genfromtxt("IDs/XopAG/XopAG_KEGG.txt", dtype = 'str')
XopAK = np.genfromtxt("IDs/XopAK/XopAK_KEGG.txt", dtype = 'str')
XopAO1 = np.genfromtxt("IDs/XopAO1/XopAO1_KEGG.txt", dtype = 'str')
XopC2 = np.genfromtxt("IDs/XopC2/XopC2_KEGG.txt", dtype = 'str')
XopE1 = np.genfromtxt("IDs/XopE1/XopE1_KEGG.txt", dtype = 'str')
XopE4 = np.genfromtxt("IDs/XopE4/XopE4_KEGG.txt", dtype = 'str')
XopK = np.genfromtxt("IDs/XopK/XopK_KEGG.txt", dtype = 'str')
XopL = np.genfromtxt("IDs/XopL/XopL_KEGG.txt", dtype = 'str')
XopN = np.genfromtxt("IDs/XopN/XopN_KEGG.txt", dtype = 'str')
XopQ = np.genfromtxt("IDs/XopQ/XopQ_KEGG.txt", dtype = 'str')
XopR = np.genfromtxt("IDs/XopR/XopR_KEGG.txt", dtype = 'str')
XopV = np.genfromtxt("IDs/XopV/XopV_KEGG.txt", dtype = 'str')
XopZ = np.genfromtxt("IDs/XopV/XopV_KEGG.txt", dtype = 'str')

#Llamado del método para obtener los primeros ec
primeros_ec_AvrBs2 = obtener_primer_ec(AvrBs2)
primeros_ec_XB3 = obtener_primer_ec(XB3)
primeros_ec_XopAE = obtener_primer_ec(XopAE)
primeros_ec_XopAG = obtener_primer_ec(XopAG)
primeros_ec_XopAK = obtener_primer_ec(XopAK)
primeros_ec_XopAO1 = obtener_primer_ec(XopAO1)
primeros_ec_XopC2 = obtener_primer_ec(XopC2)
primeros_ec_XopE1 = obtener_primer_ec(XopE1)
primeros_ec_XopE4 = obtener_primer_ec(XopE4)
primeros_ec_XopK = obtener_primer_ec(XopK)
primeros_ec_XopL = obtener_primer_ec(XopL)
primeros_ec_XopN = obtener_primer_ec(XopN)
primeros_ec_XopQ = obtener_primer_ec(XopQ)
primeros_ec_XopR = obtener_primer_ec(XopR)
primeros_ec_XopV = obtener_primer_ec(XopV)
primeros_ec_XopZ = obtener_primer_ec(XopZ)

#Llamado del método para obtener los tipos de enzimas
tipos_enzimas_avrbs2 = asignar_tipo_enzima(primeros_ec_AvrBs2)
tipos_enzimas_XB3 = asignar_tipo_enzima(primeros_ec_XB3)
tipos_enzimas_XopAE = asignar_tipo_enzima(primeros_ec_XopAE)
tipos_enzimas_XopAG = asignar_tipo_enzima(primeros_ec_XopAG)
tipos_enzimas_XopAK = asignar_tipo_enzima(primeros_ec_XopAK)
tipos_enzimas_XopAO1 = asignar_tipo_enzima(primeros_ec_XopAO1)
tipos_enzimas_XopC2 = asignar_tipo_enzima(primeros_ec_XopC2)
tipos_enzimas_XopE1 = asignar_tipo_enzima(primeros_ec_XopE1)
tipos_enzimas_XopE4 = asignar_tipo_enzima(primeros_ec_XopE4)
tipos_enzimas_XopK = asignar_tipo_enzima(primeros_ec_XopK)
tipos_enzimas_XopL = asignar_tipo_enzima(primeros_ec_XopL)
tipos_enzimas_XopN = asignar_tipo_enzima(primeros_ec_XopN)
tipos_enzimas_XopQ = asignar_tipo_enzima(primeros_ec_XopQ)
tipos_enzimas_XopR = asignar_tipo_enzima(primeros_ec_XopR)
tipos_enzimas_XopV = asignar_tipo_enzima(primeros_ec_XopV)
tipos_enzimas_XopZ = asignar_tipo_enzima(primeros_ec_XopZ)  

#Llamado del código para conteo de enzimas
df_conteo_avrbs2 = contar_tipos_enzimas(primeros_ec_AvrBs2)
df_conteo_XB3 = contar_tipos_enzimas(primeros_ec_XB3)
df_conteo_XopAE = contar_tipos_enzimas(primeros_ec_XopAE)
df_conteo_XopAG = contar_tipos_enzimas(primeros_ec_XopAG)
df_conteo_XopAK = contar_tipos_enzimas(primeros_ec_XopAK)
df_conteo_XopAO1 = contar_tipos_enzimas(primeros_ec_XopAO1)
df_conteo_XopC2 = contar_tipos_enzimas(primeros_ec_XopC2)
df_conteo_XopE1 = contar_tipos_enzimas(primeros_ec_XopE1)
df_conteo_XopE4 = contar_tipos_enzimas(primeros_ec_XopE4)
df_conteo_XopK = contar_tipos_enzimas(primeros_ec_XopK)
df_conteo_XopL = contar_tipos_enzimas(primeros_ec_XopL)
df_conteo_XopN = contar_tipos_enzimas(primeros_ec_XopN)
df_conteo_XopQ = contar_tipos_enzimas(primeros_ec_XopQ)
df_conteo_XopR = contar_tipos_enzimas(primeros_ec_XopR)
df_conteo_XopV = contar_tipos_enzimas(primeros_ec_XopV)
df_conteo_XopZ = contar_tipos_enzimas(primeros_ec_XopZ)

#Creación de gráficas de tipos de enzimas para cada efector
sns.set(style="darkgrid")
datos_avrbs2 = pd.value_counts(tipos_enzimas_avrbs2)
datos_avrbs2.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: AvrBs2', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/AvrBs2/AvrBs2_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight')   

sns.set(style="darkgrid")
datos_XB3 = pd.value_counts(tipos_enzimas_XB3)
datos_XB3.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XB3', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XB3/XB3_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight')   

sns.set(style="darkgrid")
datos_XopAE = pd.value_counts(tipos_enzimas_XopAE)
datos_XopAE.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopAE', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopAE/XopAE_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight')   

sns.set(style="darkgrid")
datos_XopAG = pd.value_counts(tipos_enzimas_XopAG)
datos_XopAG.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopAG', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopAG/XopAG_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight')  

sns.set(style="darkgrid")
datos_XopAK = pd.value_counts(tipos_enzimas_XopAK)
datos_XopAK.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopAK', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopAK/XopAK_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight')  

sns.set(style="darkgrid")
datos_XopAO1 = pd.value_counts(tipos_enzimas_XopAO1)
datos_XopAO1.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopAO1', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopAO1/XopAO1_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight')  

sns.set(style="darkgrid")
datos_XopC2 = pd.value_counts(tipos_enzimas_XopC2)
datos_XopC2.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopC2', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopC2/XopC2_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight')  

sns.set(style="darkgrid")
datos_XopE1 = pd.value_counts(tipos_enzimas_XopE1)
datos_XopE1.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopE1', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopE1/XopE1_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight')  

sns.set(style="darkgrid")
datos_XopE4 = pd.value_counts(tipos_enzimas_XopE4)
datos_XopE4.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopE4', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopE4/XopE4_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight') 

sns.set(style="darkgrid")
datos_XopK = pd.value_counts(tipos_enzimas_XopK)
datos_XopK.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopK', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopK/XopK_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight') 

sns.set(style="darkgrid")
datos_XopL = pd.value_counts(tipos_enzimas_XopL)
datos_XopL.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopL', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopL/XopL_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight') 

sns.set(style="darkgrid")
datos_XopN = pd.value_counts(tipos_enzimas_XopN)
datos_XopN.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopN', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopN/XopN_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight') 

sns.set(style="darkgrid")
datos_XopQ = pd.value_counts(tipos_enzimas_XopQ)
datos_XopQ.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopQ', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopQ/XopQ_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight') 

sns.set(style="darkgrid")
datos_XopR = pd.value_counts(tipos_enzimas_XopR)
datos_XopR.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopR', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopR/XopR_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight') 

sns.set(style="darkgrid")
datos_XopV = pd.value_counts(tipos_enzimas_XopV)
datos_XopV.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopV', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopV/XopV_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight') 

sns.set(style="darkgrid")
datos_XopZ = pd.value_counts(tipos_enzimas_XopV)
datos_XopZ.sort_index().plot.barh()
tfont = {'fontname':'Times New Roman'}
plt.title('Tipos de enzimas: XopZ', **tfont)
plt.xlabel('Tipos enzimáticos', **tfont)
plt.ylabel('Cantidad', **tfont)
pl.savefig('IDs/XopZ/XopZ_tipos_enzimas.png', format='png', dpi=900, bbox_inches='tight') 

kok = pd.read_csv("IDs/AvrBs2/query_AvrBs2.ko", sep = '\t', header = None)
interactor_raw_intento = pd.read_csv("IDs/AvrBs2/AvrBs2_ID_RAW.txt",sep = '\t', header = None)
interactor_raw_intento
no_encontrados = np.genfromtxt("IDs/AvrBs2/no_encontrado_AvrBs2.txt", dtype = 'str')
intento=crear_data_frame_manid_ko_araid_ec(interactor_raw_intento, kok, no_encontrados)
print(intento)
len(intento)
print(len(intento)==len(interactor_raw_intento))
type(kok[1][25])

#Llamado de método para llenar cada archivo de query con 0´s por cada NaN
cambiar_query("IDs/AvrBs2/query_AvrBs2.ko")
cambiar_query("IDs/XB3/query_XB3.ko")
cambiar_query("IDs/XopAE/query_XopAE.ko")
cambiar_query("IDs/XopAG/query_XopAG.ko")
cambiar_query("IDs/XopAK/query_XopAK.ko")
cambiar_query("IDs/XopAO1/query_XopAO1.ko")
cambiar_query("IDs/XopC2/query_XopC2.ko")
cambiar_query("IDs/XopE1/query_XopE1.ko")
cambiar_query("IDs/XopE4/query_XopE4.ko")
cambiar_query("IDs/XopK/query_XopK.ko")
cambiar_query("IDs/XopL/query_XopL.ko")
cambiar_query("IDs/XopN/query_XopN.ko")
cambiar_query("IDs/XopQ/query_XopQ.ko")
cambiar_query("IDs/XopR/query_XopR.ko")
cambiar_query("IDs/XopV/query_XopV.ko")
cambiar_query("IDs/XopZ/query_XopZ.ko")



#importe de la matriz de enzimas


df_cluster_funciones = pd.read_csv("cluster_funciones.csv",  sep = ';', index_col = 0)
df_cluster_rutas = pd.read_csv("cluster_rutas.csv", sep =  ';', index_col = 0)

#heatmap#
index = []
for i in df_cluster_rutas.index:
    index.append(i)
print(index)

hm_cluster_funciones = sns.heatmap(df_cluster_funciones, annot = False, square = False, cmap = 'gist_heat')
plt.title("Enzymes by cellular function")
plt.savefig('hm_cluster_funciones.pdf', dpi = 1000, bbox_inches = 'tight')
      

sns.set(font_scale = 0.55)
hm_cluster_rutas = sns.heatmap(df_cluster_rutas, annot = False, yticklabels = index, cmap = 'gist_heat')
plt.title("Enzymes by KAAS pathway type")
plt.savefig('hm_cluster_rutas.pdf', dpi = 1000, bbox_inches = 'tight')

#Importar base de datos cassavacyc#
df_cassava_cyc = pd.read_csv("bases_manihot/cassvacyc_plantcyc_pathways_manihot_esculenta.txt", sep = "\t")

df_raw = pd.read_csv("/home/gustavo/Documents/proyecto_tesis/Datos_crudos/DESeq2/analisis_de_expresión_diferencial2_sin_filtro_de_por_baja_expresion_(33032 genes)/DESeq2_Diff_expression_all_comparisons_All genes without filter.csv")
    
    
df_raw.head() 
df_raw.columns.values    

df_raw.rename(columns = {'Row-names': 'id_manihot'}, inplace = True)
prefijo = 'Manes.'
nuevo_id = []
for x in df_raw['id_manihot']:
    x = prefijo + x[6:len(x)]
    nuevo_id.append(x)

df_raw['id_manihot'] = nuevo_id


#Llamado del método para crear los archivos con id_manihot, id_arabidopsis, ec y ko, algunos interactores serán repetidos
df_ko_AvrBs2 = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/AvrBs2/AvrBs2_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/AvrBs2/query_AvrBs2.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/AvrBs2/no_encontrado_AvrBs2.txt", dtype = 'str'))
df_ko_AvrBs2.head()

df_ko_XB3 = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XB3/XB3_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XB3/query_XB3.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XB3/no_encontrado_XB3.txt", dtype = 'str'))
df_ko_XB3.head()

df_ko_XopAE = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopAE/XopAE_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopAE/query_XopAE.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XopAE/no_encontrado_XopAE.txt", dtype = 'str'))
df_ko_XopAE.head()

df_ko_XopAG = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopAG/XopAG_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopAG/query_XopAG.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XopAG/no_encontrado_XopAG.txt", dtype = 'str'))
df_ko_XopAG.head()

df_ko_XopAK = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopAK/XopAK_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopAK/query_XopAK.ko", sep = '\t', header = None),
                                                   np.genfromtxt("IDs/XopAK/no_encontrado_XopAK.txt", dtype = 'str'))
df_ko_XopAK.head()

df_ko_XopAO1 = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopAO1/XopAO1_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopAO1/query_XopAO1.ko", sep = '\t', header = None),
                                                  np.loadtxt("IDs/XopAO1/no_encontrado_XopAO1.txt", dtype = 'str'))
df_ko_XopAO1.head()

df_ko_XopC2 = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopC2/XopC2_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopC2/query_XopC2.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XopC2/no_encontrado_XopC2.txt", dtype = 'str'))
df_ko_XopC2.head()

df_ko_XopE1 = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopE1/XopE1_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopE1/query_XopE1.ko", sep = '\t', header = None),
                                                   np.genfromtxt("IDs/XopE1/no_encontrado_XopE1.txt", dtype = 'str'))
df_ko_XopE1.head()

df_ko_XopE4 = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopE4/XopE4_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopE4/query_XopE4.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XopE4/no_encontrado_XopE4.txt", dtype = 'str'))
df_ko_XopE4.head()

df_ko_XopK = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopK/XopK_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopK/query_XopK.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XopK/no_encontrado_XopK.txt", dtype = 'str'))
df_ko_XopK.head()

df_ko_XopL = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopL/XopL_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopL/query_XopL.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XopN/no_encontrado_XopN.txt", dtype = 'str'))
df_ko_XopL.head()

df_ko_XopN = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopN/XopN_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopN/query_XopN.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XopN/no_encontrado_XopN.txt", dtype = 'str'))
df_ko_XopN.head()

df_ko_XopQ = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopQ/XopQ_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopQ/query_XopQ.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XopQ/no_encontrado_XopQ.txt", dtype = 'str'))
df_ko_XopQ.head()

df_ko_XopR = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopR/XopR_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopR/query_XopR.ko", sep = '\t', header = None),
                                                   np.genfromtxt("IDs/XopR/no_encontrado_XopR.txt", dtype = 'str'))
df_ko_XopR.head()

df_ko_XopV = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopV/XopV_ID_RAW.txt", sep = '\t',header = None), 
                                                  pd.read_csv("IDs/XopV/query_XopV.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XopV/no_encontrado_XopV.txt", dtype = 'str'))
df_ko_XopV.head()

df_ko_XopZ = crear_data_frame_manid_ko_araid_ec(pd.read_csv("IDs/XopZ/XopZ_ID_RAW.txt", sep = '\t',header = None),
                                                  pd.read_csv("IDs/XopZ/query_XopZ.ko", sep = '\t', header = None),
                                                  np.genfromtxt("IDs/XopZ/no_encontrado_XopZ.txt", dtype = 'str'))
df_ko_XopZ.head()

lista_dfs_guardar = [df_ko_AvrBs2, df_ko_XB3, df_ko_XopAE, df_ko_XopAG, df_ko_XopAK, df_ko_XopAO1, df_ko_XopC2, 
                     df_ko_XopE1, df_ko_XopE4, df_ko_XopK, df_ko_XopL, df_ko_XopN, df_ko_XopQ, df_ko_XopR, df_ko_XopV,
                     df_ko_XopZ]

#añadir columnas de log2foldchange y valor p por medio de un crossref en la base original y cambiar nombre de las columnas 
df_ko_AvrBs2 = pd.merge(df_ko_AvrBs2, df_raw[['id_manihot', 'AvrBs2-Control___log2FoldChange','AvrBs2-Control___padj']],
                       how = 'left', on = 'id_manihot')
df_ko_AvrBs2.rename(columns = {'AvrBs2-Control___log2FoldChange':'log2_fold_change',
                               'AvrBs2-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XB3 = pd.merge(df_ko_XB3, df_raw[['id_manihot', 'XB3-Control___log2FoldChange', 'XB3-Control___padj']], 
                     on = 'id_manihot', how = 'left')
df_ko_XB3.rename(columns = {'XB3-Control___log2FoldChange': 'log2_fold_change',
                            'XB3-Control___padj':'valor_p_ajustado'}, inplace = True)
df_ko_XopAE = pd.merge(df_ko_XopAE, df_raw[['id_manihot', 'XopAE-Control___log2FoldChange', 'XopAE-Control___padj']],
                       on = 'id_manihot', how = 'left')
df_ko_XopAE.rename(columns = { 'XopAE-Control___log2FoldChange': 'log2_fold_change',
                              'XopAE-Control___padj':'valor_p_ajustado'}, inplace = True)
df_ko_XopAG = pd.merge(df_ko_XopAG, df_raw[['id_manihot', 'XopAG-Control___log2FoldChange', 'XopAG-Control___padj']],
                       on = 'id_manihot', how = 'left')
df_ko_XopAG.rename(columns = { 'XopAG-Control___log2FoldChange': 'log2_fold_change',
                              'XopAG-Control___padj':'valor_p_ajustado'}, inplace = True)
df_ko_XopAK = pd.merge(df_ko_XopAK, df_raw[['id_manihot', 'XopAK-Control___log2FoldChange', 'XopAK-Control___padj']],
                       on = 'id_manihot', how = 'left')
df_ko_XopAK.rename(columns = {'XopAK-Control___log2FoldChange':'log2_fold_change',
                               'XopAK-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopAO1.columns.values
df_ko_XopAO1 = pd.merge(df_ko_XopAO1, df_raw[['id_manihot', 'XopAO1-Control___log2FoldChange', 'XopAO1-Control___padj']],
                        on = 'id_manihot', how = 'left')
df_ko_XopAO1.rename(columns = {'XopAO1-Control___log2FoldChange':'log2_fold_change',
                               'XopAO1-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopC2 = pd.merge(df_ko_XopC2, df_raw[['id_manihot', 'XopC2-Control___log2FoldChange', 'XopC2-Control___padj']],
                       on = 'id_manihot', how = 'left')
df_ko_XopC2.rename(columns = {'XopC2-Control___log2FoldChange':'log2_fold_change',
                               'XopC2-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopE1 =pd.merge(df_ko_XopE1, df_raw[['id_manihot', 'XopE1-Control___log2FoldChange', 'XopE1-Control___padj']],
                      on = 'id_manihot', how = 'left')
df_ko_XopE1.rename(columns = {'XopE1-Control___log2FoldChange':'log2_fold_change',
                               'XopE1-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopE4 = pd.merge(df_ko_XopE4, df_raw[['id_manihot', 'XopE4-Control___log2FoldChange', 'XopE4-Control___padj']],
                       on = 'id_manihot', how = 'left')
df_ko_XopE4.rename(columns = {'XopE4-Control___log2FoldChange':'log2_fold_change',
                               'XopE4-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopK = pd.merge(df_ko_XopAK, df_raw[['id_manihot', 'XopK-Control___log2FoldChange', 'XopK-Control___padj']],
                      on = 'id_manihot', how = 'left')
df_ko_XopK.rename(columns = {'XopK-Control___log2FoldChange':'log2_fold_change',
                               'XopK-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopL = pd.merge(df_ko_XopL, df_raw[['id_manihot', 'XopL-Control___log2FoldChange', 'XopL-Control___padj']],
                      on = 'id_manihot', how = 'left')
df_ko_XopL.rename(columns = {'XopL-Control___log2FoldChange':'log2_fold_change',
                               'XopL-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopN = pd.merge(df_ko_XopN, df_raw[['id_manihot', 'XopN-Control___log2FoldChange',  'XopN-Control___padj']],
                      on = 'id_manihot', how = 'left')
df_ko_XopN.rename(columns = {'XopN-Control___log2FoldChange':'log2_fold_change',
                               'XopN-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopQ = pd.merge(df_ko_XopQ, df_raw[['id_manihot', 'XopQ-Control___log2FoldChange', 'XopQ-Control___padj']],
                      on = 'id_manihot', how = 'left')
df_ko_XopQ.rename(columns = {'XopQ-Control___log2FoldChange':'log2_fold_change',
                               'XopQ-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopR = pd.merge(df_ko_XopR, df_raw[['id_manihot', 'XopR-Control___log2FoldChange', 'XopR-Control___padj']],
                      on = 'id_manihot', how = 'left')
df_ko_XopR.rename(columns = {'XopR-Control___log2FoldChange':'log2_fold_change',
                               'XopR-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopV = pd.merge(df_ko_XopV, df_raw[['id_manihot','XopV-Control___log2FoldChange',   'XopV-Control___padj']],
                      on = 'id_manihot', how = 'inner')
df_ko_XopV.rename(columns = {'XopV-Control___log2FoldChange':'log2_fold_change',
                               'XopV-Control___padj': 'valor_p_ajustado'}, inplace = True)
df_ko_XopZ = pd.merge(df_ko_XopZ, df_raw[['id_manihot', 'XopZ-Control___log2FoldChange',  'XopZ-Control___padj']],
                      on = 'id_manihot', how = 'inner')
df_ko_XopZ.rename(columns = {'XopZ-Control___log2FoldChange':'log2_fold_change',
                               'XopZ-Control___padj': 'valor_p_ajustado'}, inplace = True)

#Alistar los ids de cada efector como una sola lista#
lista_bds = [df_ko_AvrBs2, df_ko_XB3, df_ko_XopAE, df_ko_XopAG, df_ko_XopAK, df_ko_XopAO1, df_ko_XopC2, 
                     df_ko_XopE1, df_ko_XopE4, df_ko_XopK, df_ko_XopL, df_ko_XopN, df_ko_XopQ, df_ko_XopR, df_ko_XopV,
                     df_ko_XopZ]
lista_bds
#hacer el crossref para cada lista en la base de casavacyc#
resultados_crossref_cassavacyc_ecs = []
for x in lista_bds:
    nuevo_crossref = crossref(x['id_manihot'], df_cassava_cyc, 'id_manihot', 'EC', unir = False)
    resultados_crossref_cassavacyc_ecs.append(nuevo_crossref)
resultados_crossref_cassavacyc_ecs
    
for x in lista_bds:
    x['EC'] = añadir_prefijo_comas(x['EC'], "EC-")

lista_resultados_comparacion_ec = []
for x, y in zip(resultados_crossref_cassavacyc_ecs, lista_bds):
    crossref_nuevo = crossref(x['id_manihot'], y, 'id_manihot', 'EC', unir = False)
    union = x.assign(EC_crossref = crossref_nuevo['EC'])
    columna_comparacion = np.where(union['EC'] == union['EC_crossref'], True, False)
    final = union.assign(comparacion = columna_comparacion)
    lista_resultados_comparacion_ec.append(final)
lista_resultados_comparacion_ec



df_ko_AvrBs2.to_csv("IDs/AvrBs2/df_ko_AvrBs2.txt", index = False)
df_ko_XB3.to_csv("IDs/XB3/df_ko_XB3.txt", index = False)
df_ko_XopAE.to_csv("IDs/XopAE/df_ko_XopAE.txt", index = False)
df_ko_XopAG.to_csv("IDs/XopAG/df_ko_XopAG.txt", index = False)
df_ko_XopAK.to_csv("IDs/XopAK/df_ko_XopAK.txt", index = False)
df_ko_XopAO1.to_csv("IDs/XopAO1/df_ko_XopAO1.txt", index = False)
df_ko_XopC2.to_csv("IDs/XopC2/df_ko_XopC2.txt", index = False)
df_ko_XopE1.to_csv("IDs/XopE1/df_ko_XopE1.txt", index = False)
df_ko_XopE4.to_csv("IDs/XopE4/df_ko_XopE4.txt", index = False)
df_ko_XopK.to_csv("IDs/XopK/df_ko_XopK.txt", index = False)
df_ko_XopL.to_csv("IDs/XopL/df_ko_XopL.txt", index = False)
df_ko_XopN.to_csv("IDs/XopN/df_ko_XopN.txt", index = False)
df_ko_XopQ.to_csv("IDs/XopQ/df_ko_XopQ.txt", index = False)
df_ko_XopR.to_csv("IDs/XopR/df_ko_XopR.txt", index = False)
df_ko_XopV.to_csv("IDs/XopV/df_ko_XopV.txt", index = False)
df_ko_XopZ.to_csv("IDs/XopZ/df_ko_XopZ.txt", index = False)


#Análisis de clusters#

df_enzymes = pd.read_csv("analisis_rutas/cluster_rutas.csv", index_col=0, sep = ";")
df_enzymes_norm = pd.DataFrame()
keys = df_enzymes.columns.values
for i in keys:
    df_enzymes_norm[i] = df_enzymes[i]/df_enzymes[i].sum()
df_enzymes_norm

"""Puntaje silueta para decidir número de clusters"""
k_posibles = []
scores = []
for k in range (2, 20):
    clusterer = KMeans(n_clusters = k)
    preds = clusterer.fit_predict(df_enzymes_norm)
    score = silhouette_score(df_enzymes_norm, preds)
    k_posibles.append(k)
    scores.append(score)
sns.set(style = "darkgrid")
plt.plot(k_posibles, scores)

"""Análisis de clusters como tal"""
k = 4
modelo_1 = KMeans(n_clusters = k)
modelo_1.fit(df_enzymes_norm)
num_cluster_1_puntos = modelo_1.labels_.tolist()
num_cluster_1_puntos


clusters_rutas = []
for i in range(0, 4):
    cluster_id = df_enzymes_norm[modelo_1.labels_ == i].index
    clusters_rutas.append(cluster_id)
clusters_rutas














