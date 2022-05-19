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
        #O retorno da função anterior poderia ser uma lista, mas isso exigiria alterações
        imageList = periodCollection.toList(periodCollection.size())
        for i in range( periodCollection.size().getInfo() ):
            image = ee.Image(imageList.get(i))
            baselineList.append(image)

    gc.collect

    #Reempacotar a lista de imagens em uma coleção
    repackCollection = ee.ImageCollection.fromImages( baselineList )

    return repackCollection

#--------------------------------------------
def timeseries_postfire(startData,\
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

    beginingYMD = datetime.strptime(startData, "%Y-%m-%d") #mais antiga (inicio do período total)
    
    #Final da série histórica
    endYMD = datetime.strptime(endData, "%Y-%m-%d") #mais recente (fim)
    startYMD = endYMD - relativedelta.relativedelta(days=lapse)

    #Período de análise -- construção da série de dados
    period_pack = []
    period_index = 1
    while True:

        #Seleção da coleção do período  >>> Fazer T=0 aqui! na primeira vez =0 vai dar problema?
        periodCollection = period_collection(startYMD,endYMD,1,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem)

        ##Cálculo de índices espectrais (do período) e do "máximo" (do período)
        imagePeriod, periodCollection = include_attributes(periodCollection, \
                                                           sensorName, \
                                                           listAtts, \
                                                           str(period_index), \
                                                           'median')#'max')

        #Incluir a imagem extraída do período em uma nova coleção (a qual será processada ao fim)
        #...juntamente com a imagem "baseline-prefogo"
        period_pack.append(imagePeriod)

        #Aplicar um coletor de lixo aqui?
        gc.collect()

        #Testar essa condição (a ordem pode estar trocada!)
        if ((startYMD - relativedelta.relativedelta(days=lapse)) - beginingYMD).days <= 0:
            break #Condição para quebrar o while...
        else:
            endYMD = startYMD
            startYMD = endYMD - relativedelta.relativedelta(days=lapse)
            period_index += 1
    
    #Reempacotamento da coleção, após ajuste de nuvens/background
    repackPeriodCollection = ee.ImageCollection.fromImages( period_pack )
    
    return repackPeriodCollection

#--------------------------------------------
def period_collection(iniData,endData,T,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem):

    info = specIndex.sensor_info(sensorName['multispec'])

    _startYDM = iniData - relativedelta.relativedelta(months=T)
    startData = str(_startYDM.year)+'-'+str(_startYDM.month)+'-'+str(_startYDM.day)
    str_endData =   str(endData.year)+'-'+str(endData.month)+'-'+str(endData.day)
    
    pre_collection = ee.ImageCollection(info['completeSensorName']) \
                       .filterBounds(geometry) \
                       .filterDate(startData,str_endData) \
                       .filterMetadata('CLOUD_COVER','less_than', cloudCoverMaxPercentage)

    pre_collection, vecTime = maskingFunctions.collection_chrono_sort(pre_collection)

    collection, LoI = maskingFunctions.repack_collection(pre_collection,       #Coleção sob análise
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
        return period_collection(iniData,endData,T+1,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem)

#--------------------------------------------
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

#--------------------------------------------
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

#--------------------------------------------
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

#--------------------------------------------
def build_base_image(prefireMedian,
               postfireSerie,
               listAtts,
               baseline_sufix,
               statType,
               alpha,
               beta,
               geometry,
               scale,
               sensorName,   #pode ser desnecessário...
               coordSystem):

    #Transformar "postFireSeries" em uma imagem com as bandas que são relevantes
    baseImage =  bild_time_stack_image(postfireSerie,listAtts)
    
    #Junção das imagens em uma única base
    for item in listAtts:
        str_includeBand = item+'_'+baseline_sufix
        baseImage = baseImage.addBands( prefireMedian.select(str_includeBand) )

    #Cálculo dos deltaNBRs via "band math do EE"
    baseDeltaImage = compute_delta_nbr(baseImage,baseline_sufix)

    #Adicionar estatísticas dos atributos no período
    baseDeltaImage = compute_stats_atts(baseImage,baseDeltaImage,listAtts,statType)

    #Operar sobre os deltaNBRs? (PRECISA DE AJUSTES!? o beta está estranho...)
    baseDeltaImage = compute_burn_att(baseImage,baseDeltaImage,alpha,beta)

    #Temp...
    print(baseDeltaImage.clipToBoundsAndScale(geometry=geometry,scale=scale).select('sum_thres_deltaNBR').getThumbUrl({'min': 0, 'max': 10}))

    return baseDeltaImage
    
#--------------------------------------------
def build_base_image__var(prefireMedian,
               postfireSerie,
               #baseImage,
               listAtts,
               baseline_sufix,
               statType,
               alpha,
               beta,
               early_fire,
               geometry,
               scale,
               sensorName,   #pode ser desnecessário...
               coordSystem):

    ##Transformar "postFireSeries" em uma imagem com as bandas que são relevantes
    baseImage =  bild_time_stack_image(postfireSerie,listAtts)
    
    #Junção das imagens em uma única base
    for item in listAtts:
        str_includeBand = item+'_'+baseline_sufix
        baseImage = baseImage.addBands( prefireMedian.select(str_includeBand) )

    #Cálculo dos deltaNBRs via "band math do EE"
    baseDeltaImage = compute_delta_nbr(baseImage,baseline_sufix)

    #Adicionar estatísticas dos atributos no período
    baseDeltaImage = compute_stats_atts(baseImage,baseDeltaImage,listAtts,statType)

    #Operar sobre os deltaNBRs? (PRECISA DE AJUSTES!? o beta está estranho...)
    baseDeltaImage = compute_burn_att(baseImage,baseDeltaImage,alpha,beta)
    
    #Incluir as bandas de referência/atributos
    if early_fire == 1:
        baseDeltaImage = compute_atts_early_fire(baseImage,baseDeltaImage,listAtts,geometry,scale,coordSystem)   
    if early_fire == 2:
        baseDeltaImage = compute_atts_early_first_fire(baseImage,baseDeltaImage,listAtts,geometry,scale,coordSystem)   
    if early_fire == 3:
        baseDeltaImage = compute_atts_early_first_fire_22nov21(baseImage,baseDeltaImage,listAtts,geometry,scale,coordSystem)   
    if early_fire == 4:
        baseDeltaImage = compute_atts_early_first_fire_22nov21_modis(baseImage,baseDeltaImage,listAtts,geometry,scale,coordSystem)   

    return baseDeltaImage

#--------------------------------------------
def build_base_image__varModis(prefireMedian,
               postfireSerie,
               listAtts,
               baseline_sufix,
               statType,
               alpha,
               beta,
               early_fire,
               geometry,
               scale,
               coordSystem,
               baseModis):

    ##Transformar "postFireSeries" em uma imagem com as bandas que são relevantes
    baseImage =  bild_time_stack_image(postfireSerie,listAtts)
    
    #Junção das imagens em uma única base
    for item in listAtts:
        str_includeBand = item+'_'+baseline_sufix
        baseImage = baseImage.addBands( prefireMedian.select(str_includeBand) )

    #Cálculo dos deltaNBRs via "band math do EE"
    baseDeltaImage = compute_delta_nbr(baseImage,baseline_sufix)

    #Adicionar estatísticas dos atributos no período
    baseDeltaImage = compute_stats_atts(baseImage,baseDeltaImage,listAtts,statType)

    #Operar sobre os deltaNBRs? (PRECISA DE AJUSTES!? o beta está estranho...)
    #baseDeltaImage = compute_burn_att(baseImage,baseDeltaImage,alpha,beta)
    baseDeltaImage = compute_burn_att_modisSupport(baseImage,baseDeltaImage,baseModis,alpha,beta)
    
    #Incluir as bandas de referência/atributos
    if early_fire == 1:
        baseDeltaImage = compute_atts_early_fire(baseImage,baseDeltaImage,listAtts,geometry,scale,coordSystem)   
    if early_fire == 2:
        baseDeltaImage = compute_atts_early_first_fire(baseImage,baseDeltaImage,listAtts,geometry,scale,coordSystem)   
    if early_fire == 3:
        baseDeltaImage = compute_atts_early_first_fire_22nov21(baseImage,baseDeltaImage,listAtts,geometry,scale,coordSystem)   
    
    return baseDeltaImage

#--------------------------------------------
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

#--------------------------------------------
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

#--------------------------------------------
def bild_time_stack_image(postfireSerie,listAtts):

    baseImage = ee.Image()

    #Para cada imagem da coleção "postfireSerie"
    dim = postfireSerie.size().getInfo()
    imageList = postfireSerie.toList(postfireSerie.size())
    for i in range(dim):
        #imageInfo = ee.Image(imageList.get(i)).getInfo() #???   #Removido em 22nov21
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

#--------------------------------------------
def bild_time_stack_cloud_masks(postfireSerie):

    strMasks = 'cloudAndBackgroundMask'
    baseMaskImage = ee.Image()

    #Para cada imagem da coleção "postfireSerie"
    dim = postfireSerie.size().getInfo()
    imageList = postfireSerie.toList(postfireSerie.size())
    for i in range(dim):
        #image = ee.Image(imageList.get(i))
        #temp = image.select(strMasks)
        temp = ee.Image(imageList.get(i)).select(strMasks).rename(strMasks+'_'+str(i+1))
        baseMaskImage = baseMaskImage.addBands(temp)

    return baseMaskImage

#--------------------------------------------
def compute_stats_atts(baseImage,baseDeltaImage,listAtts,statType):
    
    bandList = baseImage.bandNames().getInfo()
    for atts in listAtts:
        listSelBands = []
        exp = re.compile(atts+'_\d')
        
        for item in bandList:
            out = exp.match(item)
            if out != None: listSelBands.append(item)
        
        if len(listSelBands) > 0:
            if statType == 'mean':
                temp = baseImage.select(listSelBands).reduce(ee.Reducer.mean()).rename(atts+'_'+statType)
            if statType == 'median':
                temp = baseImage.select(listSelBands).reduce(ee.Reducer.median()).rename(atts+'_'+statType)
            if statType == 'stdev':
                temp = baseImage.select(listSelBands).reduce(ee.Reducer.sampleStdDev()).rename(atts+'_'+statType)
            
            baseDeltaImage = baseDeltaImage.addBands(temp)

    return baseDeltaImage

#--------------------------------------------
def compute_atts_early_fire(baseImage,baseDeltaImage,listAtts,geometry,scale,coordSystem):
    
    #Criação da lista de limiares-deltaNBR
    bandDeltaList = baseDeltaImage.bandNames().getInfo()
    listSelThres = []
    exp = re.compile('thres_deltaNBR_\d')
    for item in bandDeltaList:
            out = exp.match(item)
            if out != None: 
                listSelThres.append(item)

    #Criação da lista de atributos-instantes
    bandList = baseImage.bandNames().getInfo()
    dictAttsList = {}
    for att in listAtts:
        exp = re.compile(att+'_\d')
        tempList = []
        for item in bandList:         
            out = exp.match(item)
            if out != None: 
                tempList.append(item)
        dictAttsList[att] = tempList

    #Criação de uma lista de imagens, para cada atributo, que deve conter os valores a cada último evento de fogo
    dictBeforeImages = {}
    for att in listAtts:
        refAtt = att+'_'+'ref'
        bandNameAtt = dictAttsList[att]
        indexNameAtt = [int(i.split('_')[1]) for i in bandNameAtt]
        
        #Imagem de referência, inicialmente toda nula (apenas a base para receber os valores dos atributos em "caso de fogo")
        B = baseImage.expression('band * 0',{'band': baseImage.select(att+'_1')}).rename(refAtt)      
        
        for item in listSelThres:

            indexThres = int(item.split('_')[2])
            posAntAtt = np.where( np.array(indexNameAtt) == (indexThres -1) )[0]

            if len(posAntAtt) > 0:
                selAtt = att+'_'+str( indexNameAtt[posAntAtt[0]] )

                AuxImg = ee.Image().addBands(
                                            [B.select(refAtt),baseDeltaImage.select(item),baseImage.select(selAtt)] 
                                            ).rename(['Null','Ref','Thres','Att'])
 
                dictOper = {'Ref': AuxImg.select('Ref'), 'Thres': AuxImg.select('Thres'), 'Att':  AuxImg.select('Att')}
                expression = '( !(Thres > 0)*Ref ) + ( (Thres > 0)*Att )'
                temp = AuxImg.expression(expression,dictOper).rename('updateRef')

                B = B.addBands(temp)
                B = B.select('updateRef').rename(refAtt)
        
        dictBeforeImages[refAtt] = B
    
    #Adicionar as bandas de referencia ao dataFrame principal?
    for newRefBand in dictBeforeImages:
        baseDeltaImage = baseDeltaImage.addBands( dictBeforeImages[newRefBand] )

    return baseDeltaImage

#--------------------------------------------
def compute_atts_early_first_fire(baseImage,baseDeltaImage,listAtts,geometry,scale,coordSystem):
    
    #Criação da lista de limiares-deltaNBR
    bandDeltaList = baseDeltaImage.bandNames().getInfo()
    listSelThres = []
    exp = re.compile('thres_deltaNBR_\d')
    for item in bandDeltaList:
            out = exp.match(item)
            if out != None: 
                listSelThres.append(item)

    #Criação da lista de atributos-instantes
    bandList = baseImage.bandNames().getInfo()
    dictAttsList = {}
    for att in listAtts:
        exp = re.compile(att+'_\d')
        tempList = []
        for item in bandList:         
            out = exp.match(item)
            if out != None: 
                tempList.append(item)
        dictAttsList[att] = tempList

    #Criação de uma lista de imagens, para cada atributo, que deve conter os valores a cada último evento de fogo
    dictBeforeImages = {}
    for att in listAtts:
        refAtt = att+'_'+'ref'
        bandNameAtt = dictAttsList[att]
        indexNameAtt = [int(i.split('_')[1]) for i in bandNameAtt]
        
        #Imagem de referência, inicialmente toda nula (apenas a base para receber os valores dos atributos em "caso de fogo")
        B = baseImage.expression('band * 0 - 9999.99',{'band': baseImage.select(att+'_1')}).rename(refAtt) #-9999.99 é um dummy!
        for item in listSelThres:

            indexThres = int(item.split('_')[2])
            #posAntAtt = np.where( np.array(indexNameAtt) == (indexThres -1)[0]
            posAntAtt = np.where( np.array(indexNameAtt) == (indexThres -1) )[0]

            if len(posAntAtt) > 0:
                selAtt = att+'_'+str( indexNameAtt[posAntAtt[0]] )

                AuxImg = ee.Image().addBands(
                                            [B.select(refAtt),baseDeltaImage.select(item),baseImage.select(selAtt)] 
                                            ).rename(['Null','Ref','Thres','Att'])
 
                dictOper = {'Ref': AuxImg.select('Ref'), 'Thres': AuxImg.select('Thres'), 'Att':  AuxImg.select('Att'), 'Base':  AuxImg.select('Base')}
                #expression = '( !(Thres > 0)*Ref ) + ( (Thres > 0)*Att )'
                expression = '(Ref == -9999.99) * ( !(Thres > 0)*Ref ) + ( (Thres > 0)*Att )'
                temp = AuxImg.expression(expression,dictOper).rename('updateRef')

                B = B.addBands(temp)
                B = B.select('updateRef').rename(refAtt)
        
        dictBeforeImages[refAtt] = B
    
    #Adicionar as bandas de referencia ao dataFrame principal?
    for newRefBand in dictBeforeImages:
        baseDeltaImage = baseDeltaImage.addBands( dictBeforeImages[newRefBand] )

    return baseDeltaImage

#--------------------------------------------
def compute_atts_early_first_fire_22nov21(baseImage,baseDeltaImage,listAtts,geometry,scale,coordSystem):
    
    #Criação da lista de limiares-deltaNBR
    bandDeltaList = baseDeltaImage.bandNames().getInfo()
    listSelThres = []
    exp = re.compile('thres_deltaNBR_\d')
    for item in bandDeltaList:
            out = exp.match(item)
            if out != None: 
                listSelThres.append(item)

    #Criação da lista de atributos-instantes
    bandList = baseImage.bandNames().getInfo()
    dictAttsList = {}
    for att in listAtts:
        exp = re.compile(att+'_\d')
        tempList = []
        for item in bandList:         
            out = exp.match(item)
            if out != None: 
                tempList.append(item)
        dictAttsList[att] = tempList

    #Criação de uma lista de imagens, para cada atributo, que deve conter os valores a cada último evento de fogo
    dictBeforeImages = {}
    for att in listAtts:
        refAtt = att+'_'+'ref'
        bandNameAtt = dictAttsList[att]
        indexNameAtt = [int(i.split('_')[1]) for i in bandNameAtt]
        
        #Imagem de referência, inicialmente toda nula (apenas a base para receber os valores dos atributos em "caso de fogo")
        Support = baseImage.expression('band * 0',{'band': baseImage.select(att+'_1')}).rename(refAtt) #-9999.99 é um dummy!
        pack = []
        for item in listSelThres:
            indexThres = int(item.split('_')[2])
            posAntAtt = np.where( np.array(indexNameAtt) == (indexThres -1) )[0]

            if len(posAntAtt) > 0:
                selAtt = att+'_'+str( indexNameAtt[posAntAtt[0]] )

                AuxImg = ee.Image().addBands(
                                            [Support.select(refAtt),baseDeltaImage.select(item),baseImage.select(selAtt)] 
                                            ).rename(['Null','Ref','Thres','Att'])
 
                dictOper = {'Ref': AuxImg.select('Ref'), 'Thres': AuxImg.select('Thres'), 'Att':  AuxImg.select('Att'), 'Base':  AuxImg.select('Base')}
                expression = '((Ref <= 0) * (Thres > 0)) * Att' #Para extracao dos atributos locais
                expSupport = 'Ref + ((Ref <= 0) * (Thres > 0))'       #Para identificacao dos locais disponíveis
                temp = AuxImg.expression(expression,dictOper).rename('updateRef')
                Support = AuxImg.expression(expSupport,dictOper).rename(refAtt)

                pack.append( temp )
        
        dictBeforeImages[refAtt] = ee.ImageCollection.fromImages( pack ).sum().rename(refAtt)
    
    #Adicionar as bandas de referencia ao dataFrame principal?
    for newRefBand in dictBeforeImages:
        baseDeltaImage = baseDeltaImage.addBands( dictBeforeImages[newRefBand] )

    return baseDeltaImage

#--------------------------------------------
def compute_burn_att(baseImage,baseDeltaImage,alpha,beta):

    bandList = baseDeltaImage.bandNames().getInfo()
    listSelBands = []
    listSelBandsNDWI = []
    exp = re.compile('deltaNBR_\d')
    for item in bandList:
            out = exp.match(item)
            if out != None: 
                listSelBands.append(item)
                listSelBandsNDWI.append('NDWI_' + item.split('_')[1]) #Se existe 'deltaNBR_\d', deve exitir p NDWI

    #Cálculo dos limiares sobre cada "deltaNBR" e adição à imagem "base"
    for item, ndwi in zip(listSelBands,listSelBandsNDWI):
        dictOper = {'alpha': alpha, 'beta': beta, 'band': baseDeltaImage.select(item), 'water':  baseImage.select(ndwi)}
        #temp = baseDeltaImage.expression('band >= alpha',dictOper).rename('thres_'+item)
        temp = baseDeltaImage.expression('(band >= alpha) * (water < beta)',dictOper).rename('thres_'+item)
        baseDeltaImage = baseDeltaImage.addBands(temp)

    #Cálculo de uma imagem que soma todas as ocorrências de indicativo de queimada (acima de alpha)
    bandList2 = baseDeltaImage.bandNames().getInfo()
    listSelThresBands = []
    exp = re.compile('thres_deltaNBR_\d')
    for item in bandList2:
            out = exp.match(item)
            if out != None: listSelThresBands.append(item)

    if len(listSelThresBands) > 0:
        temp = baseDeltaImage.select(listSelThresBands).reduce(ee.Reducer.sum()).rename('sum_thres_deltaNBR')
        baseDeltaImage = baseDeltaImage.addBands(temp) #temp contém uma imagem com a posição dos locais atingidos por fogo

    return baseDeltaImage

#--------------------------------------------
def compute_burn_att_modisSupport(baseImage,baseDeltaImage,baseModis,alpha,beta):

    bandList = baseDeltaImage.bandNames().getInfo()
    listSelBands = []
    listSelBandsNDWI = []
    exp = re.compile('deltaNBR_\d')
    for item in bandList:
            out = exp.match(item)
            if out != None: 
                listSelBands.append(item)
                listSelBandsNDWI.append('NDWI_' + item.split('_')[1]) #Se existe 'deltaNBR_\d', deve exitir p NDWI

    #Cálculo dos limiares sobre cada "deltaNBR" e adição à imagem "base"    
    for item, ndwi in zip(listSelBands,listSelBandsNDWI):
        dictOper = {'alpha': alpha, 'beta': beta, \
                    'band': baseDeltaImage.select(item), 'water':  baseImage.select(ndwi), \
                    'modis': baseModis.select('focos')}
        
        temp = baseDeltaImage.expression('(band >= alpha) * (water < beta) * modis',dictOper).rename('thres_'+item)
        baseDeltaImage = baseDeltaImage.addBands(temp)

    #Cálculo de uma imagem que soma todas as ocorrências de indicativo de queimada (acima de alpha)
    bandList2 = baseDeltaImage.bandNames().getInfo()
    listSelThresBands = []
    exp = re.compile('thres_deltaNBR_\d')
    for item in bandList2:
            out = exp.match(item)
            if out != None: listSelThresBands.append(item)

    if len(listSelThresBands) > 0:
        temp = baseDeltaImage.select(listSelThresBands).reduce(ee.Reducer.sum()).rename('sum_thres_deltaNBR')
        baseDeltaImage = baseDeltaImage.addBands(temp) #temp contém uma imagem com a posição dos locais atingidos por fogo

    return baseDeltaImage

#--------------------------------------------
def build_base_dataframe(baseDeltaImage,
                         listAtts,
                         statType,
                         geometry,
                         scale,
                         coordSystem):

    #Iniciar construção do DataFrame
    d = {}
    defaultDummy = 0.0 #-9999.0
    base = ee.Image(defaultDummy).blend( baseDeltaImage )

    lat, lon, value = ext_lat_lon_pixel(base, geometry, ['sum_thres_deltaNBR'], scale, coordSystem)
    d['Latitude'] = lat
    d['Longitude'] = lon
    d['sum_thres_deltaNBR'] = value[0]

    #Inclusão dos atributos do período
    listaRefs = [item+'_'+statType for item in listAtts]
    _, _, values = ext_lat_lon_pixel(base, geometry, listaRefs, scale, coordSystem)
    for item, index in zip(listaRefs,range(len(listaRefs))):
        d[item] = values[index]
        
    gc.collect() #Coletor de lixo 

    #Salvar a imagem para conferência
    tab = pd.DataFrame()
    tab = tab.from_dict(d)

    return tab

#--------------------------------------------
def build_base_dataframe_UFD(baseImage,
                             listAtts,
                             #statType,
                             geometry,
                             scale,
                             coordSystem):

    #Iniciar construção do DataFrame
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
        nbrW = transformModule.nbr_weight(value[ind])
        d['W_'+att] = nbrW
        bandsTiff.append('W_'+att)

    #Cálculo dos atributos remapeados
    for singleAtt in listAtts:
        listSelAtt, listSelIndexAtt = transformModule.get_attribute_index_from_list(singleAtt,bandNames)

        #Fase de mapeamento dos atributos
        if singleAtt == 'NDVI':
            for name, item in zip(listSelAtt,listSelIndexAtt):
                d['map_'+name] = value[item]
                bandsTiff.append( 'map_'+name )

        if singleAtt == 'NBR':
            for name, item in zip(listSelAtt,listSelIndexAtt):
                d['map_'+name] = value[item]
                bandsTiff.append( 'map_'+name )

    
    gc.collect() #Coletor de lixo 

    #Construir e preencher o dataFrame resultante
    tab = pd.DataFrame()
    tab = tab.from_dict(d)

    #Salvar a imagem para conferência
    import saveModule
    path_out = '/home/rogerio/Desktop/saidaUFD_SI.tif'
    saveModule.save_tiff_from_df(tab,bandsTiff,0.0,path_out,coordSystem)

    return tab


#--------------------------------------------
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

#--------------------------------------------
def early_assessment_period(startData,\
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

    #mais antiga (inicio do período total)
    beginingYMD = datetime.strptime(startData, "%Y-%m-%d") 
    
    #Final da série histórica
    endYMD = datetime.strptime(endData, "%Y-%m-%d") #mais recente (fim)
    startYMD = endYMD - relativedelta.relativedelta(days=lapse)

    #Período de análise -- construção da série de dados
    period_pack = []
    period_index = 1
    while True:

        #Seleção da coleção do período
        TT=0
        periodCollection = period_collection(startYMD,endYMD,TT,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem)

        ##Cálculo de índices espectrais (do período) e do "máximo" (do período)
        imagePeriod, periodCollection = include_attributes(periodCollection, \
                                                           sensorName, \
                                                           listAtts, \
                                                           str(period_index), \
                                                           'mean')

        
        period_pack.append(imagePeriod)
        gc.collect()

        if ((startYMD - relativedelta.relativedelta(days=lapse)) - beginingYMD).days <= 0:
            break #Condição para quebrar o while...
        else:
            endYMD = startYMD
            startYMD = endYMD - relativedelta.relativedelta(days=lapse)
            period_index += 1
    
    #Reempacotamento da coleção, após ajuste de nuvens/background
    repackPeriodCollection = ee.ImageCollection.fromImages( period_pack )
    
    return repackPeriodCollection

#--------------------------------------------
def build_assessment_image(assessmentSerie,
                           listAtts,
                           statType):
                           

    #Transformar "postFireSeries" em uma imagem com as bandas que são relevantes
    baseImage =  bild_time_stack_image(assessmentSerie,listAtts)
    
    #Adicionar estatísticas dos atributos no período
    baseImage = compute_stats_atts(baseImage,baseImage,listAtts,statType)

    
    return baseImage

#--------------------------------------------
def build_assessment_dataframe(baseDeltaImage,
                               listAtts,
                               statType,
                               geometry,
                               scale,
                               coordSystem):

    #Iniciar construção do DataFrame
    d = {}
    defaultDummy = 0.0 #-9999.0
    #base = ee.Image(defaultDummy).blend(ee.Image( baseDeltaImage ))
    base = ee.Image(defaultDummy).blend( baseDeltaImage )

    #Inclusão dos atributos do período
    for item in listAtts:
        item += '_'+statType
        lat, lon, value = ext_lat_lon_pixel(base, geometry, [item], scale, coordSystem)
        d[item] = value[0]

    d['Latitude'] = lat
    d['Longitude'] = lon

    gc.collect() #Coletor de lixo 

    #Construir e preencher o dataFrame resultante
    tab = pd.DataFrame()
    tab = tab.from_dict(d)

    return tab

#--------------------------------------------
def build_assessment_early_instants_image(assessmentSerie,
                           listAtts):
                           

    #Transformar "postFireSeries" em uma imagem com as bandas que são relevantes
    baseImage =  bild_time_stack_image(assessmentSerie,listAtts)

    #Obter as informações de nuvens e background em cada imagem da série de avaliação
    cloudMasksInstants =  bild_time_stack_cloud_masks(assessmentSerie)

    instantsAtts = baseImage.bandNames().getInfo()
    
    t = []
    for item in instantsAtts:
        if listAtts[0] in item: t.append(item)
    refIndex = [np.int64(i.split('_')[1]) for i in t]

    dictImageInstants = {} #Dicionário com as imagens de atributos em cada instante de "Assessment Period"

    #Para cada instante...
    for ind in refIndex:
        tempList = []
        for att in listAtts:
            tempList.append( att+'_'+str(ind) )
        dictImageInstants[ind] = ee.Image().addBands( baseImage.select(tempList) )

    return baseImage, dictImageInstants, cloudMasksInstants

#--------------------------------------------
def build_assessment_instant_dataframe(baseDeltaImage,
                                      instant,
                                      listAtts,
                                      geometry,
                                      scale,
                                      coordSystem):


    #Iniciar construção do DataFrame
    d = {}
    defaultDummy = 0.0 #-9999.0
    base = ee.Image(defaultDummy).blend( baseDeltaImage )

    #Inclusão dos atributos do período
    for item in listAtts:
        lat, lon, value = ext_lat_lon_pixel(base, geometry, [item+'_'+str(instant)], scale, coordSystem)
        d[item] = value[0]

    d['Latitude'] = lat
    d['Longitude'] = lon

    gc.collect() #Coletor de lixo 

    #Construir e preencher o dataFrame resultante
    tab = pd.DataFrame()
    tab = tab.from_dict(d)

    return tab

#--------------------------------------------
def build_assessment_instant_dataframe__withClouds(baseDeltaImage,
                                      cloudMasksInstants,
                                      instant,
                                      listAtts,
                                      geometry,
                                      scale,
                                      coordSystem):


    #Iniciar construção do DataFrame
    d = {}
    defaultDummy = 0.0 #-9999.0
    base = ee.Image(defaultDummy).blend( baseDeltaImage )

    cloudBand = cloudMasksInstants.select('cloudAndBackgroundMask_'+str(instant)).rename('cloud')
    base = base.addBands(cloudBand)

    #---------------------------------------otimizando
    list = ['cloud']
    listName = ['cloud']
    for item in listAtts:
        list.append( item+'_'+str(instant) )
        listName.append( item )
    
    lat, lon, values = ext_lat_lon_pixel(base, geometry, list, scale, coordSystem)    
    d['Latitude'] = lat
    d['Longitude'] = lon
    
    #for item,index in zip(list,range(len(list))):
    for item,index in zip(listName,range(len(listName))):
        d[item] = values[index]

    gc.collect() #Coletor de lixo

    #Construir e preencher o dataFrame resultante
    tab = pd.DataFrame()
    tab = tab.from_dict(d)

    return tab

#--------------------------------------------
def build_assessment_instant_dataframe__withClouds_andRefModis(baseDeltaImage,
                                      cloudMasksInstants,
                                      instant,
                                      listAtts,
                                      geometry,
                                      scale,
                                      coordSystem,refModis):


    #Iniciar construção do DataFrame
    d = {}
    defaultDummy = 0.0 
    base = ee.Image(defaultDummy).blend( baseDeltaImage )

    cloudBand = cloudMasksInstants.select('cloudAndBackgroundMask_'+str(instant)).rename('cloud')
    base = base.addBands(cloudBand)
    refModisBand = refModis.select(['focos']).rename('refModis')
    base = base.addBands(refModisBand)

    #---------------------------------------otimizando
    list = ['cloud','refModis']
    listName = ['cloud','refModis']
    for item in listAtts:
        list.append( item+'_'+str(instant) )
        listName.append( item )
    
    #Formar uma lista geral com os atributos 
    lat, lon, values = ext_lat_lon_pixel(base, geometry, list, scale, coordSystem)    

    #Iniciar dicionário
    d['Latitude'] = lat
    d['Longitude'] = lon
    
    for item,index in zip(listName,range(len(listName))):
        d[item] = values[index]

    gc.collect() #Coletor de lixo 

    #Construir e preencher o dataFrame resultante
    tab = pd.DataFrame()
    tab = tab.from_dict(d)

    return tab