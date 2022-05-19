
#--------------------------------------------
#IMPORTS ------------------------------------
from datetime import datetime
from dateutil import relativedelta
import gc 
import numpy as np
import pandas as pd
import re

#--------------------------------------------
#GEE ----------------------------------------
import ee
ee.Initialize() 

#--------------------------------------------
#PERSONAL MODULES ---------------------------
import specIndex
import maskingFunctions
import transformModule
import saveModule
import moduleUtils


#==================================================================
#RACIONAL:
def baseline_prefire(startData,\
                     endData,\
                     geometry,\
                     scale,\
                     spamMedianSupport,\
                     UsefulnessThres,\
                     sensorName,\
                     cloudCoverMaxPercentage,\
                     coordSystem):

    #Extrai os dados para um dado período específico startData~endData
    startYMD = datetime.strptime(startData, "%Y-%m-%d") #mais antigo (inicio)
    endYMD = datetime.strptime(endData, "%Y-%m-%d")     #mais recente (fim)
    relative = relativedelta.relativedelta(endYMD,startYMD)
    deltaMonths = 12*relative.years + relative.months

    baselineList = []
    for delta in range(0,deltaMonths+1):
        refIni = startYMD + relativedelta.relativedelta(months=delta)
        refDataIni = str(refIni.year)+'-'+str(refIni.month)+'-'+str(refIni.day)
        print("Searching in "+refDataIni+" reference date...")

        #Desempacotar a série para correção de nuvens e reempacotar na sequência
        periodCollection = expand_period_collection(refDataIni,1,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem)
        
        #Agrupar o resultado em uma coleção ainda maior, que vai compor todo o "baseline period"
        imageList = periodCollection.toList(periodCollection.size())
        for i in range( periodCollection.size().getInfo() ):
            image = ee.Image(imageList.get(i))
            baselineList.append(image)

    gc.collect

    #Reempacotar a lista de imagens em uma coleção
    repackCollection = ee.ImageCollection.fromImages( baselineList )

    return repackCollection



#==================================================================
def baseline_prefire_norec(startData,\
                     endData,\
                     geometry,\
                     scale,\
                     spamMedianSupport,\
                     UsefulnessThres,\
                     sensorName,\
                     cloudCoverMaxPercentage,\
                     coordSystem):


    boolFlag, periodCollection, LoI = expand_period_collection_norec(startData,endData,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem)

    baselineList = []
    if boolFlag:
        #Agrupar o resultado em uma coleção ainda maior, que vai compor todo o "baseline period"
        imageList = periodCollection.toList(periodCollection.size())
        for i in range( periodCollection.size().getInfo() ):
            image = ee.Image(imageList.get(i))
            baselineList.append(image)
        
        repackCollection = ee.ImageCollection.fromImages( baselineList )
        gc.collect
        return True, repackCollection, LoI
    else:
        print('Problema... não encontrou imagem boa no periodo')
        return False, 'Null', []      




#==================================================================
def timeseries_postfire_simple(startData,\
                        endData,\
                        #lapse,\
                        listAtts,\
                        geometry,\
                        scale,\
                        spamMedianSupport,\
                        UsefulnessThres,\
                        sensorName,\
                        cloudCoverMaxPercentage,\
                        coordSystem):

    #Período de análise -- construção da série de dados
    startYMD = datetime.strptime(startData, "%Y-%m-%d") #mais recente (fim) #Final da série histórica
    endYMD = datetime.strptime(endData, "%Y-%m-%d") #mais recente (fim)

    #Seleção da coleção do período
    flag, periodCollection, LoI = period_collection_norec(startYMD,endYMD,0,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem)

    if flag:
        #Calcular os atributos para cada uma das imagens da coleção e definir um índice para estas!
        period_pack = []
        instant_list = []
        dim = periodCollection.size().getInfo()
        imageList = periodCollection.toList(dim)
        for ind in range(dim):
            image = ee.Image(imageList.get(ind))

            XXX = ee.ImageCollection.fromImages( [image] )
            imagePeriod, periodCollection = include_attributes_nostats(XXX, \
                                                               sensorName, \
                                                               listAtts, \
                                                               str(ind+1), \
                                                               'median')
            period_pack.append(imagePeriod)
            instant =  datetime.strptime(image.getInfo()['properties']['SENSING_TIME'].split('T')[0], "%Y-%m-%d")
            #instant_list.append( instant )
            instant_list.append( LoI[ind] )

            #Aplicar um coletor de lixo aqui?
            gc.collect()
    
        #Reempacotamento da coleção, após ajuste de nuvens/background
        repackPeriodCollection = ee.ImageCollection.fromImages( period_pack )

        return flag, repackPeriodCollection, instant_list
    else:
        return False, 'nodata', 'nodata', []



#==================================================================
def timeseries_postfire_simple__25set21(startData,\
                        endData,\
                        lapse,\
                        listAtts,\
                        geometry,\
                        scale,\
                        spamMedianSupport,\
                        UsefulnessThres,\
                        sensorName,\
                        cloudCoverMaxPercentage,\
                        coordSystem):

    #Período de análise -- construção da série de dados
    startYMD = datetime.strptime(startData, "%Y-%m-%d") #mais recente (fim) #Final da série histórica
    endYMD = datetime.strptime(endData, "%Y-%m-%d") #mais recente (fim)

    #Seleção da coleção do período
    flag, periodCollection, vecChrono = period_collection_norec(startYMD,endYMD,0,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem)

    if flag:
        #Calcular os atributos para cada uma das imagens da coleção e definir um índice para estas!
        period_pack = []
        instant_list = []
        dim = periodCollection.size().getInfo()
        imageList = periodCollection.toList(dim)
        for ind in range(dim):
            image = ee.Image(imageList.get(ind))

            XXX = ee.ImageCollection.fromImages( [image] )
            imagePeriod, periodCollection = include_attributes_nostats(XXX, \
                                                               sensorName, \
                                                               listAtts, \
                                                               str(ind+1), \
                                                               'median')
            period_pack.append(imagePeriod)
            instant =  datetime.strptime(image.getInfo()['properties']['SENSING_TIME'].split('T')[0], "%Y-%m-%d")
            instant_list.append( instant )

            #Aplicar um coletor de lixo aqui?
            gc.collect()
    
        #Reempacotamento da coleção, após ajuste de nuvens/background
        repackPeriodCollection = ee.ImageCollection.fromImages( period_pack )

        return flag, instant_list, repackPeriodCollection, vecChrono
    else:
        return False, 'nodata', 'nodata', 'nodata'


#==================================================================
def period_collection_norec(iniData,endData,T,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem):

    info = specIndex.sensor_info(sensorName['multispec'])

    _startYDM = iniData - relativedelta.relativedelta(months=T)
    startData = str(_startYDM.year)+'-'+str(_startYDM.month)+'-'+str(_startYDM.day)
    str_endData =   str(endData.year)+'-'+str(endData.month)+'-'+str(endData.day)
    
    pre_collection = ee.ImageCollection(info['completeSensorName']) \
                       .filterBounds(geometry) \
                       .filterDate(startData,str_endData) \
                       .filterMetadata('CLOUD_COVER','less_than', cloudCoverMaxPercentage)

    pre_collection, init_LoI = maskingFunctions.collection_chrono_sort(pre_collection)

    collection, LoI = maskingFunctions.repack_collection(pre_collection,            #Coleção sob análise
                                                    geometry,                  #Geometria da coleção
                                                    scale,                     #Escala dos dados
                                                    spamMedianSupport,         #Período (meses) usados na geração de uma base de suporte (img. mediana)
                                                    sensorName,                #Nome do sensor (útil a criação das máscaras de nuvem)
                                                    cloudCoverMaxPercentage,   #Percentual de cobertura de nuvem
                                                    UsefulnessThres)           #Utilidade da imagem (>50%)

    found = collection.size().getInfo()

    #Chamada recursiva
    #Enquanto não encontra imagens válidas no intervalo (refData-T ~ refData)
    #o valor de T aumenta mais, expandindo o intervalo de busca para (refData-(T+1) ~ refData)
    if found == 0:
        return False, 'null', []
    else:
        #Chamada recursiva...
        #return True, period_collection(iniData,endData,T+1,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem)
        return True, collection, LoI



#==================================================================
def period_collection_norec__25set21(iniData,endData,T,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem):

    info = specIndex.sensor_info(sensorName['multispec'])

    _startYDM = iniData - relativedelta.relativedelta(months=T)
    startData = str(_startYDM.year)+'-'+str(_startYDM.month)+'-'+str(_startYDM.day)
    str_endData =   str(endData.year)+'-'+str(endData.month)+'-'+str(endData.day)
    
    pre_collection = ee.ImageCollection(info['completeSensorName']) \
                       .filterBounds(geometry) \
                       .filterDate(startData,str_endData) \
                       .filterMetadata('CLOUD_COVER','less_than', cloudCoverMaxPercentage)

    pre_collection, vecChrono = maskingFunctions.collection_chrono_sort(pre_collection)

    collection = maskingFunctions.repack_collection(pre_collection,            #Coleção sob análise
                                                    geometry,                  #Geometria da coleção
                                                    scale,                     #Escala dos dados
                                                    spamMedianSupport,         #Período (meses) usados na geração de uma base de suporte (img. mediana)
                                                    sensorName,                #Nome do sensor (útil a criação das máscaras de nuvem)
                                                    cloudCoverMaxPercentage,   #Percentual de cobertura de nuvem
                                                    UsefulnessThres)           #Utilidade da imagem (>50%)

    found = collection.size().getInfo()

    #Chamada recursiva
    #Enquanto não encontra imagens válidas no intervalo (refData-T ~ refData)
    #o valor de T aumenta mais, expandindo o intervalo de busca para (refData-(T+1) ~ refData)
    if found == 0:
        return False, 'null', 'null'
    else:
        #Chamada recursiva...
        #return True, period_collection(iniData,endData,T+1,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem)
        return True, collection, vecChrono        



#==================================================================
def include_attributes(collection, sensorName, listAtts, sufix, stat):

    #Cálculo/inclusão de índices espectrais
    for item in listAtts:
        if item == 'NDVI':
            collection = specIndex.add_ndvi_collection(collection,sensorName['multispec'])
        if item == 'NDWI':
            #collection = specIndex.add_ndwi_collection(collection,sensorName['multispec'])
            collection = specIndex.add_mndwi_collection(collection,sensorName['multispec'])
        if item == 'SAVI':
            collection = specIndex.add_savi_collection(collection,sensorName['multispec'])
        if item == 'NBR':
            collection = specIndex.add_nbr_collection(collection,sensorName['multispec'])

    #Cálculo da estatística da coleção
    if stat == 'median': statImage = collection.median()
    if stat == 'mean': statImage = collection.mean()
    if stat == 'max': statImage = collection.max()
    if stat == 'min': statImage = collection.min()    

    #Renomear as bandas após a estatísticas
    statImage = rename_bands(statImage,sufix,listAtts)

    return statImage, collection


#==================================================================
def include_attributes_nostats(collection, sensorName, listAtts, sufix, stat):

    #Cálculo/inclusão de índices espectrais
    for item in listAtts:
        if item == 'NDVI':
            collection = specIndex.add_ndvi_collection(collection,sensorName['multispec'])
        if item == 'NDWI':
            #collection = specIndex.add_ndwi_collection(collection,sensorName['multispec'])
            collection = specIndex.add_mndwi_collection(collection,sensorName['multispec'])
        if item == 'SAVI':
            collection = specIndex.add_savi_collection(collection,sensorName['multispec'])
        if item == 'NBR':
            collection = specIndex.add_nbr_collection(collection,sensorName['multispec'])

    dim = collection.size().getInfo()
    if dim > 1:
        print('...')

    #Cálculo da estatística da coleção
    if stat == 'median': statImage = collection.median()
    if stat == 'mean': statImage = collection.mean()
    if stat == 'max': statImage = collection.max()
    if stat == 'min': statImage = collection.min()    

    #Renomear as bandas após a estatísticas
    statImage = rename_bands(statImage,sufix,listAtts)

    return statImage, collection


#==================================================================
def rename_bands(statImage,sufix,listAtts):

    bandList = statImage.bandNames().getInfo()
    newList = []
    for band in bandList:

        if (band.lower() in listAtts) or (band.upper() in listAtts):
            newList.append( band.upper()+'_'+sufix )
        else:
            newList.append( band ) 
    
    statImage = statImage.select(bandList).rename(newList)

    return statImage


#==================================================================
def build_base_dataframe_UFD(baseImage,
                             listAtts,
                             #statType,
                             geometry,
                             scale,
                             coordSystem,path_out,LoI):

    #...iniciar construção do DataFrame
    d = {}
    bandsTiff = []

    defaultDummy = 0.0 #-9999.0
    base = ee.Image(defaultDummy).blend( baseImage )

    bandNames = base.bandNames().getInfo()

    #Usar Regex para mapear os atributos de interesse/potencial para o mapeamento
    listPostAtts = transformModule.select_attributes_from_list(listAtts+['deltaNBR'],bandNames)

    #Formar uma lista geral com os atributos de "listPostDeltaNBR" e "listPostAtts"
    lat, lon, value = ext_lat_lon_pixel(base, geometry, listPostAtts, scale, coordSystem)

    #Início do dicionário...
    d['Latitude'] = lat
    d['Longitude'] = lon

    listSelAttDeltaNBR, listSelIndexDeltaNBR = transformModule.get_attribute_index_from_list('deltaNBR',listPostAtts)

    #Cálculo do peso a partir do deltaNBR
    for att,ind in zip(listSelAttDeltaNBR,listSelIndexDeltaNBR):
        #nbrW = transformModule.nbr_weight(value[ind])
        #nbrW = transformModule.delta_nbr_weight(value[ind])
        #d['W_'+att] = nbrW
        d['W_'+att] = value[ind]
        bandsTiff.append('W_'+att)


    #Cálculo dos atributos remapeados
    for singleAtt in listAtts:
        #listSelAtt, listSelIndexAtt = transformModule.get_attribute_index_from_list(singleAtt,bandNames)
        listSelAtt, listSelIndexAtt = transformModule.get_attribute_index_from_list(singleAtt,listPostAtts)

        #Fase de mapeamento dos atributos
        if singleAtt == 'NDVI':
            for name, item in zip(listSelAtt,listSelIndexAtt):
                #d['map_'+name] = transformModule.ndvi_map(value[item])
                d['map_'+name] = value[item]
                bandsTiff.append( 'map_'+name )

        if singleAtt == 'NBR':
            for name, item in zip(listSelAtt,listSelIndexAtt):
                #d['map_'+name] = transformModule.nbr_map(value[item])
                d['map_'+name] = value[item]
                bandsTiff.append( 'map_'+name )

    
    gc.collect() #Coletor de lixo :: Precisa desalocar a var. "base" ou outra qualquer?

    #Construct and fill in the resulting dataFrame
    tab = pd.DataFrame()
    tab = tab.from_dict(d)

    #Salvar a imagem para conferir...
    #path_out = '/home/rogerio/Desktop/UFD/___saidaUFD_SIv11_mai-out20_1monthSpan.tif'
    saveModule.save_tiff_from_df(tab,bandsTiff,0.0,path_out,coordSystem)
    moduleUtils.save_image_names(path_out,bandsTiff)#,vecChrono)
    moduleUtils.add_instant_info(path_out,LoI)
    
    return tab



#==================================================================
def build_base_dataframe_UFD__25set21(baseImage,
                             listAtts,
                             #statType,
                             geometry,
                             scale,
                             coordSystem,
                             path_out,vecChrono):

    #...iniciar construção do DataFrame
    d = {}
    bandsTiff = []

    defaultDummy = 0.0 #-9999.0
    base = ee.Image(defaultDummy).blend( baseImage )

    bandNames = base.bandNames().getInfo()

    #Usar Regex para mapear os atributos de interesse/potencial para o mapeamento
    listPostAtts = transformModule.select_attributes_from_list(listAtts+['deltaNBR'],bandNames)

    #Formar uma lista geral com os atributos de "listPostDeltaNBR" e "listPostAtts"
    lat, lon, value = ext_lat_lon_pixel(base, geometry, listPostAtts, scale, coordSystem)

    #Início do dicionário...
    d['Latitude'] = lat
    d['Longitude'] = lon

    listSelAttDeltaNBR, listSelIndexDeltaNBR = transformModule.get_attribute_index_from_list('deltaNBR',listPostAtts)

    #Cálculo do peso a partir do deltaNBR
    for att,ind in zip(listSelAttDeltaNBR,listSelIndexDeltaNBR):
        #nbrW = transformModule.nbr_weight(value[ind])
        #nbrW = transformModule.delta_nbr_weight(value[ind])
        #d['W_'+att] = nbrW
        d['W_'+att] = value[ind]
        bandsTiff.append('W_'+att)


    #Cálculo dos atributos remapeados
    for singleAtt in listAtts:
        listSelAtt, listSelIndexAtt = transformModule.get_attribute_index_from_list(singleAtt,bandNames)

        #Fase de mapeamento dos atributos
        if singleAtt == 'NDVI':
            for name, item in zip(listSelAtt,listSelIndexAtt):
                #d['map_'+name] = transformModule.ndvi_map(value[item])
                d['map_'+name] = value[item]
                bandsTiff.append( 'map_'+name )

        if singleAtt == 'NBR':
            for name, item in zip(listSelAtt,listSelIndexAtt):
                #d['map_'+name] = transformModule.nbr_map(value[item])
                d['map_'+name] = value[item]
                bandsTiff.append( 'map_'+name )

    
    gc.collect() #Coletor de lixo :: Precisa desalocar a var. "base" ou outra qualquer?

    #Construct and fill in the resulting dataFrame
    tab = pd.DataFrame()
    tab = tab.from_dict(d)

    #Salvar a imagem para conferir...
    #path_out = '/home/rogerio/Desktop/UFD/saidaUFD_SIv11_mai-out20_1monthSpan.tif'
    saveModule.save_tiff_from_df(tab,bandsTiff,0.0,path_out,coordSystem)
    moduleUtils.save_image_names(path_out,bandsTiff,vecChrono)

    return tab



#==============================================
#Adotar um esquema recursivo (que retorna a data cada vez mais, até que encontre alguma imagem)
#Contempla apenas indices espectrais (teste inicial...)
def expand_period_collection(refData,T,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem):
    
    print("Busca de nivel "+str(T))

    info = specIndex.sensor_info(sensorName['multispec'])

    ymd = datetime.strptime(refData, "%Y-%m-%d")    
    startYDM = ymd - relativedelta.relativedelta(months=T)
    startData = str(startYDM.year)+'-'+str(startYDM.month)+'-'+str(startYDM.day)
    endData = refData
    
    pre_collection = ee.ImageCollection(info['completeSensorName']) \
                       .filterBounds(geometry) \
                       .filterDate(startData,endData) \
                       .filterMetadata('CLOUD_COVER','less_than', cloudCoverMaxPercentage)

    pre_collection = maskingFunctions.collection_chrono_sort(pre_collection)

    collection = maskingFunctions.repack_collection(pre_collection,            #Coleção sob análise
                                                    geometry,                  #Geometria da coleção
                                                    scale,                     #Escala dos dados
                                                    spamMedianSupport,         #Período (meses) usados na geração de uma base de suporte (img. mediana)
                                                    sensorName,                #Nome do sensor (útil a criação das máscaras de nuvem)
                                                    cloudCoverMaxPercentage,   #Percentual de cobertura de nuvem
                                                    UsefulnessThres)           #Utilidade da imagem (>50%)

    found = collection.size().getInfo()

    #Chamada recursiva
    #Enquanto não encontra imagens válidas no intervalo (refData-T ~ refData)
    #o valor de T aumenta mais, expandindo o intervalo de busca para (refData-(T+1) ~ refData)
    if found != 0:
        return collection
    else:
        #Chamada recursiva...
        return expand_period_collection(refData,T+1,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem)



#==============================================
#Adotar um esquema recursivo (que retorna a data cada vez mais, até que encontre alguma imagem)
#Contempla apenas indices espectrais (teste inicial...)
def expand_period_collection_norec(startData,endData,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem):
    
    #refData,T,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem
    #startData,endData,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem

    print("Busca sem recursão...")

    info = specIndex.sensor_info(sensorName['multispec'])
    
    pre_collection = ee.ImageCollection(info['completeSensorName']) \
                       .filterBounds(geometry) \
                       .filterDate(startData,endData) \
                       .filterMetadata('CLOUD_COVER','less_than', cloudCoverMaxPercentage)

    pre_collection, LoI_init = maskingFunctions.collection_chrono_sort(pre_collection)

    collection, LoI = maskingFunctions.repack_collection(pre_collection,            #Coleção sob análise
                                                    geometry,                  #Geometria da coleção
                                                    scale,                     #Escala dos dados
                                                    spamMedianSupport,         #Período (meses) usados na geração de uma base de suporte (img. mediana)
                                                    sensorName,                #Nome do sensor (útil a criação das máscaras de nuvem)
                                                    cloudCoverMaxPercentage,   #Percentual de cobertura de nuvem
                                                    UsefulnessThres)           #Utilidade da imagem (>50%)

    found = collection.size().getInfo()

    #Retorno do esquema sem recursao
    if found != 0:
        return True, collection, LoI
    else:
        return False, 'Null', []




#==================================================================
def build_base_image_events(prefireMedian,
                            postfireSerie,
                            listAtts,
                            baseline_sufix):
                            
    
    ##Transformar "postFireSeries" em uma imagem com as bandas que são relevantes
    baseImage =  bild_time_stack_image(postfireSerie,listAtts)
    
    #Junção das imagens em uma única base
    for item in listAtts:
        str_includeBand = item+'_'+baseline_sufix
        baseImage = baseImage.addBands( prefireMedian.select(str_includeBand) )

    #Cálculo dos deltaNBRs via "band math do EE"
    baseDeltaImage = compute_delta_nbr(baseImage,baseline_sufix)
    
    print('...') 
    baseImage = baseImage.addBands( baseDeltaImage )

    return baseImage 



#==================================================================
def compute_delta_nbr(baseImage,baseline_sufix):
    
    #Instanciação da nova imagem
    baseDeltaImage = ee.Image()

    #Definição da banda "pre-fogo"
    baselineNBR = 'NBR_'+baseline_sufix
    
    #Lista de bandas para análise (todas)
    bandList = baseImage.bandNames().getInfo()
    
    #Definindo uma expressão regular para encontrar bandas de nome "NBR_[0~9]"
    exp = re.compile('NBR_\d')
    listPostNBR = []
    for item in bandList:
        out = exp.match(item)
        if out != None: listPostNBR.append(item)
    
    #Cálculo das bandas "deltaNBR" e adição à imagem "base" (que contém toda a informação útil do processo)
    for item in listPostNBR:
        dictOper = {'preFire' : baseImage.select(baselineNBR), 'postFire': baseImage.select(item) }
        temp = baseImage.expression('preFire - postFire',dictOper).rename('delta'+item)
        baseDeltaImage = baseDeltaImage.addBands(temp)

    return baseDeltaImage


#==================================================================
def bild_time_stack_image(postfireSerie,listAtts):

    baseImage = ee.Image()

    #Para cada imagem da coleção "postfireSerie"
    dim = postfireSerie.size().getInfo()
    imageList = postfireSerie.toList(postfireSerie.size())
    for i in range(dim):
        imageInfo = ee.Image(imageList.get(i)).getInfo() #???
        image = ee.Image(imageList.get(i))
        bandList = image.bandNames().getInfo()

        listSelBands = []
        for atts in listAtts:
            exp = re.compile(atts+'_\d') #Definindo uma expressão regular para encontrar bandas de nome "item_[0~9]"
            
            for item in bandList:
                out = exp.match(item)
                if out != None: listSelBands.append(item)
        
        if len(listSelBands) > 0:
            baseImage = baseImage.addBands( image.select(listSelBands) )

    return baseImage



#==============================================
def ext_lat_lon_pixel(image, geometria, bandas, escala, coordSystem):
    image = image.addBands(ee.Image.pixelLonLat())
    dictImg = image.select(['longitude', 'latitude']+bandas) \
                   .reduceRegion(reducer=ee.Reducer.toList(),\
                                 geometry=geometria,scale=escala,\
                                 bestEffort=True,crs=coordSystem)
    bandas_valores = []
    for banda in bandas:
        bandas_valores.append(np.array(ee.List(dictImg.get(banda)).getInfo()).astype(float))

    return np.array(ee.List(dictImg.get('latitude')).getInfo()).astype(float), np.array(ee.List(dictImg.get('longitude')).getInfo()).astype(float), bandas_valores

