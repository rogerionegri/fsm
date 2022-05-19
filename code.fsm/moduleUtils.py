
#--------------------------------------------
#IMPORTS ------------------------------------
import pandas as pd
import gc

#--------------------------------------------
#GEE ----------------------------------------
import ee
ee.Initialize() 

#--------------------------------------------
#PERSONAL MODULES ---------------------------
import saveModule
import UFD_DataBase


#==================================================================
#RACIONAL:
def save_image(baseImage,
               #listAtts,
               path_out,
               geometry,
               scale,
               coordSystem):

    #...iniciar construção do DataFrame
    d = {}

    defaultDummy = 0.0 #-9999.0
    base = ee.Image(defaultDummy).blend( baseImage )

    bandNames = base.bandNames().getInfo()

    #Formar uma lista geral com os atributos de "listPostDeltaNBR" e "listPostAtts"
    lat, lon, value = UFD_DataBase.ext_lat_lon_pixel(base, geometry, bandNames, scale, coordSystem)

    #Início do dicionário...
    d['Latitude'] = lat
    d['Longitude'] = lon

    #Cálculo dos atributos remapeados
    for singleAtt, ind in zip(bandNames,range(len(bandNames))):
        d[singleAtt] = value[ind]
    
    gc.collect() #Coletor de lixo 

    #Construct and fill in the resulting dataFrame
    tab = pd.DataFrame()
    tab = tab.from_dict(d)

    #Salvar a imagem para conferir...
    saveModule.save_tiff_from_df(tab,bandNames,0.0,path_out,coordSystem)
    save_image_names(path_out,bandNames)

    return tab

#==================================================================
def save_image_names(path_out,bandNames):
    #Arquivo de texto com o nome das bandas
    path_txt = path_out+'_bands.txt'
    file = open(path_txt,"w")
    for b in bandNames: file.write(b+'\n')

    file.close()
    return True

#==================================================================
def add_instant_info(path_out,LoI):

    #Arquivo de texto com o nome das bandas
    path_txt = path_out+'_bands.txt'
    file = open(path_txt,"a")
    file.write('\n')
    for info in LoI: file.write(info+'\n')
    file.close()

    return True

#==================================================================
def add_report_info(path_out,title,vecStr):
    
    file = open(path_out,"a")
    file.write('\n')
    file.write(title+'\n')
    for info in vecStr: file.write(info+'\n')
    file.write('\n')
    file.close()

    return True