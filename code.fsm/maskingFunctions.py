#--------------------------------------------
#IMPORTS ------------------------------------
import numpy as np
import gc

from datetime import datetime
from dateutil import relativedelta

#--------------------------------------------
#GEE ----------------------------------------
import ee
ee.Initialize() 

#--------------------------------------------
#PERSONAL MODULES ---------------------------
import specIndex

#--------------------------------------------
#Função para geração da máscara de nuvem
def build_cloud_mask(image,sensorName):
    infoMSpec = specIndex.sensor_info(sensorName['multispec'])
    
    #Landasat
    if sensorName['multispec'] == 'landsat':
        bandQA = image.select(infoMSpec['quality'])
                
        shadow = bandQA.bitwiseAnd(8).eq(0)
        cloud = bandQA.bitwiseAnd(32).eq(0)
        cloudMask = shadow.multiply(cloud)

        #---------------outra abordagem...
        #"Cloud Confidence".
        bitStartCloudConfidence = 5
        bitEndCloudConfidence = 6
        qaBitsCloudConfidence = extractQABits(bandQA, bitStartCloudConfidence, bitEndCloudConfidence)
        testCloudConfidence = qaBitsCloudConfidence.gte(2)

        #"Cloud Shadow Confidence".
        bitStartShadowConfidence = 7
        bitEndShadowConfidence = 8
        qaBitsShadowConfidence = extractQABits(bandQA, bitStartShadowConfidence, bitEndShadowConfidence)
        testShadowConfidence = qaBitsShadowConfidence.gte(2)

        maskComposite = testCloudConfidence.Or(testShadowConfidence).Not()
    
    return cloudMask

#--------------------------------------------
#Função para geração da máscara de nuvem
def build_cloud_mask_V2(image,sensorName,sensorString,geometry,scale):
        
    #Landasat
    if sensorString == 'multispec':
        infoMSpec = specIndex.sensor_info(sensorName['multispec'])
        bandQA = image.select(infoMSpec['quality'])
               
        shadow = bandQA.bitwiseAnd(8).eq(0)
        cloud = bandQA.bitwiseAnd(32).eq(0)
        cloudMask = shadow.multiply(cloud)

        #---------------outra abordagem...
        #"Cloud Confidence".
        bitStartCloudConfidence = 5
        bitEndCloudConfidence = 6
        qaBitsCloudConfidence = extractQABits(bandQA, bitStartCloudConfidence, bitEndCloudConfidence)
        testCloudConfidence = qaBitsCloudConfidence.gte(2)

        #"Cloud Shadow Confidence".
        bitStartShadowConfidence = 7
        bitEndShadowConfidence = 8
        qaBitsShadowConfidence = extractQABits(bandQA, bitStartShadowConfidence, bitEndShadowConfidence)
        testShadowConfidence = qaBitsShadowConfidence.gte(2)

        #Calculate a composite mask and apply it to the image.   
        maskComposite = testCloudConfidence.Or(testShadowConfidence).Not()

        return cloudMask
    
    #Sentinel
    if sensorString == 'multispecHR':
        infoMSpec = specIndex.sensor_info(sensorName['multispecHR'])
        tempSel = image.select(infoMSpec['SCL'])
        dictOper = {'slc': tempSel}
        strOper = '(slc == 3)||(slc == 8)||(slc == 9)||(slc == 10)||(slc == 1)||(slc == -1)||(slc == 0)'
        cloudShadowMask = tempSel.expression(strOper,dictOper).Not().rename('cloudShadowMask')
        return cloudShadowMask  

#--------------------------------------------
#Função para extração de bits da banda pixel_qa
def extractQABits(qaBand, bitStart, bitEnd):
    RADIX = 2
    numBits = bitEnd - bitStart + 1
    qaBits = qaBand.rightShift(bitStart).mod(RADIX**numBits)   #Math.pow(RADIX, numBits))
    return qaBits

#--------------------------------------------
#Função para análise e adequação da coleção
def repack_collection(collection,                #Coleção sob análise
                      geometry,                  #Geometria da coleção
                      scale,                     #Escala dos dados
                      pastMonthsMedian,          #Período (meses) usados na geração de uma base de suporte (img. mediana)
                      sensorName,                #Nome do sensor (útil a criação das máscaras de nuvem)
                      cloudCoverMaxPercentage,   #Percentual de cobertura de nuvem
                      UsefulnessThres):          #Limiar de aproveitamento da imagem (cobertura de nuvem+background)
    
    kernelRadius = 2

    #Obter informações do sensor
    info = specIndex.sensor_info(sensorName['multispec'])
    BoI = info['BoI']

    repack = []
    listOfInstants = []

    dim = collection.size().getInfo()
    imageList = collection.toList(collection.size())
    for i in range(dim):
        imageInfo = ee.Image(imageList.get(i)).getInfo()
        image = ee.Image(imageList.get(i))

        #Definir a máscara desta imagem (melhorar a parametrização)
        cloudMask = build_cloud_mask(image,sensorName)

        #Definição da máscara de nuvem+background (i.e., dados faltante)
        cloudBackgroundMask_out = ee.Image(0).blend(cloudMask).rename('cloudAndBackgroundMask')

        #Com base no comportamento de "cloudBackgroundMask_out", é possível ter uma noção do aproveitamento da imagem!
        #Obtém e transforma a máscara e uma lista 0-1 que indica a pertinencia ou não em área útil
        valuesMask = np.array( ee.List(cloudBackgroundMask_out \
                                 .reduceRegion(reducer=ee.Reducer.toList(),geometry=geometry,scale=10*scale,bestEffort=True)) \
                                 .getInfo()['cloudAndBackgroundMask'] )
        #Realiza o cálculo do percentual de elementos "1" (i.e., em área útil) em relação ao total
        Usefulness = np.count_nonzero(valuesMask == 1) / valuesMask.shape[0]
        
        if Usefulness >= UsefulnessThres:

            #------------------------------------------------
            #Definir uma mediana em relação ao instante de "imageList.get(i))
            imagingTime = imageInfo['properties']['SENSING_TIME'].split('T')[0]
            #imagingTime = imageInfo['properties']['DATE_ACQUIRED']
            endYMD = datetime.strptime(imagingTime, "%Y-%m-%d")     #mais recente (fim)      
            
            factorPast, earlyCount = 1, 0
            while earlyCount == 0:
                ref = endYMD - relativedelta.relativedelta(months=factorPast*pastMonthsMedian)  #3 é um teste...
                backTime = str(ref.year)+'-'+str(ref.month)+'-'+str(ref.day)

                #Definir a busca de modo mais conveniente aqui (melhorar a parametrização)
                earlyCollection = ee.ImageCollection(info['completeSensorName']) \
                                        .filterBounds(geometry) \
                                        .filterDate(backTime,imagingTime) \
                                        .filterMetadata('CLOUD_COVER','less_than',cloudCoverMaxPercentage)
                
                earlyCount = earlyCollection.size().getInfo()
                factorPast += 1
            #-----------------------------------------------

            #earlyCollection = collection_chrono_sort(earlyCollection) #<<<Desnecessário reordenar...

            medianImage = earlyCollection.median()
            #------------------------------------------------

            #Aplicar um filtro de dilatação (circular de raio 1)
            cloudBackgroundMask_out = cloudBackgroundMask_out.focal_max(kernel=ee.Kernel.circle(radius=kernelRadius), iterations=1)

            #Adicao da mascara de dados faltantes como banda da imagem (para fins de referencia futura)
            image = image.addBands([cloudBackgroundMask_out])

            #Aplicação das máscaras In/Out sobre as BoI's
            newBandNames = []
            for band in BoI:

                #Imagem/banda original -- aplicação da máscara para remoção dos dados ausentes/contaminados
                __bandMasked_in  = image.select(band).updateMask(cloudBackgroundMask_out).rename(band+'_mask'+'In')
            
                #Low-level median -- aplicação da máscara reversa para aproveitamento dos dados
                __bandMedianMasked_out = medianImage.select(band).updateMask(cloudBackgroundMask_out.Not()).rename(band+'median_mask'+'Out')

                #Preenchimento de dados faltantes
                __blendBand = __bandMedianMasked_out.blend(__bandMasked_in).rename(band+'_corr')

                image = image.addBands([__blendBand])

            #Renomear as bandas corrigidas e originais
            bandNames = image.bandNames().getInfo()
            #for b in bandNames:
            for b in BoI:
                posOrig = bandNames.index(b)
                posCorr = bandNames.index(b+'_corr')
                bandNames[posCorr] = bandNames[posOrig]
                bandNames[posOrig] = bandNames[posOrig]+'_original' 
            image = image.select(image.bandNames().getInfo()).rename(bandNames)

            #Armazenando a imagem para instanciar a coleção em seguida
            repack.append(image)

            #Armazenar o instante da imagem extraída da coleção
            listOfInstants.append( imagingTime )

            #Aplicar um coletor de lixo aqui?
            gc.collect()

    #Reempacotamento da coleção, após ajuste de nuvens/background
    repackCollection = ee.ImageCollection.fromImages( repack )
    
    return repackCollection, listOfInstants

#--------------------------------------------
#Função para análise e adequação da coleção
def repack_collection__25set21(collection,       #Coleção sob análise
                      geometry,                  #Geometria da coleção
                      scale,                     #Escala dos dados
                      pastMonthsMedian,          #Período (meses) usados na geração de uma base de suporte (img. mediana)
                      sensorName,                #Nome do sensor (útil a criação das máscaras de nuvem)
                      cloudCoverMaxPercentage,   #Percentual de cobertura de nuvem
                      UsefulnessThres):          #Limiar de aproveitamento da imagem (cobertura de nuvem+background)

    #Obter informações do sensor
    info = specIndex.sensor_info(sensorName['multispec'])
    BoI = info['BoI']

    repack = []
    dim = collection.size().getInfo()
    imageList = collection.toList(collection.size())
    for i in range(dim):
        imageInfo = ee.Image(imageList.get(i)).getInfo()
        image = ee.Image(imageList.get(i))

        #Definir a máscara desta imagem (melhorar a parametrização)
        cloudMask = build_cloud_mask(image,sensorName)

        #Definição da máscara de nuvem+background (i.e., dados faltante)
        cloudBackgroundMask_out = ee.Image(0).blend(cloudMask).rename('cloudAndBackgroundMask')

        #Com base no comportamento de "cloudBackgroundMask_out", é possível ter uma noção do aproveitamento da imagem!
        #Obtém e transforma a máscara e uma lista 0-1 que indica a pertinencia ou não em área útil
        valuesMask = np.array( ee.List(cloudBackgroundMask_out \
                                 .reduceRegion(reducer=ee.Reducer.toList(),geometry=geometry,scale=10*scale,bestEffort=True)) \
                                 .getInfo()['cloudAndBackgroundMask'] )
        #Realiza o cálculo do percentual de elementos "1" (i.e., em área útil) em relação ao total
        Usefulness = np.count_nonzero(valuesMask == 1) / valuesMask.shape[0]
        
        if Usefulness >= UsefulnessThres:

            #------------------------------------------------
            #Definir uma mediana em relação ao instante de "imageList.get(i))
            imagingTime = imageInfo['properties']['SENSING_TIME'].split('T')[0]
            #imagingTime = imageInfo['properties']['DATE_ACQUIRED']
            endYMD = datetime.strptime(imagingTime, "%Y-%m-%d")     #mais recente (fim)      
            ref = endYMD - relativedelta.relativedelta(months=pastMonthsMedian)  #3 é um teste...
            backTime = str(ref.year)+'-'+str(ref.month)+'-'+str(ref.day)

            #Definir a busca de modo mais conveniente aqui (melhorar a parametrização)
            earlyCollection = ee.ImageCollection(info['completeSensorName']) \
                                .filterBounds(geometry) \
                                .filterDate(backTime,imagingTime) \
                                .filterMetadata('CLOUD_COVER','less_than',cloudCoverMaxPercentage)

            #earlyCollection = collection_chrono_sort(earlyCollection) #<<<Desnecessário reordenar...

            medianImage = earlyCollection.median()
            #------------------------------------------------

            #Aplicar um filtro de dilatação (circular de raio 1)
            cloudBackgroundMask_out = cloudBackgroundMask_out.focal_max(kernel=ee.Kernel.circle(radius=1), iterations=1)

            #Adicao da mascara de dados faltantes como banda da imagem (para fins de referencia futura)
            image = image.addBands([cloudBackgroundMask_out])

            #Aplicação das máscaras In/Out sobre as BoI's
            newBandNames = []
            for band in BoI:

                #Imagem/banda original -- aplicação da máscara para remoção dos dados ausentes/contaminados
                __bandMasked_in  = image.select(band).updateMask(cloudBackgroundMask_out).rename(band+'_mask'+'In')
            
                #Low-level median -- aplicação da máscara reversa para aproveitamento dos dados
                __bandMedianMasked_out = medianImage.select(band).updateMask(cloudBackgroundMask_out.Not()).rename(band+'median_mask'+'Out')

                #Preenchimento de dados faltantes
                __blendBand = __bandMedianMasked_out.blend(__bandMasked_in).rename(band+'_corr')

                image = image.addBands([__blendBand])

            #Renomear as bandas corrigidas e originais
            try:
                bandNames = image.bandNames().getInfo()
                #for b in bandNames:
                for b in BoI:
                    posOrig = bandNames.index(b)
                    posCorr = bandNames.index(b+'_corr')
                    bandNames[posCorr] = bandNames[posOrig]
                    bandNames[posOrig] = bandNames[posOrig]+'_original' 
                    image = image.select(image.bandNames().getInfo()).rename(bandNames)

                #Armazenando a imagem para instanciar a coleção em seguida
                repack.append(image)

                #Aplicar um coletor de lixo aqui?
                gc.collect()
            except:
                print('exceção...')

    #Reempacotamento da coleção, após ajuste de nuvens/background
    repackCollection = ee.ImageCollection.fromImages( repack )
    
    return repackCollection

#--------------------------------------------
#Função para análise e adequação da coleção
def repack_collection_sentinel(collection,                #Coleção sob análise
                               geometry,                  #Geometria da coleção
                               scale,                     #Escala dos dados
                               pastMonthsMedian,          #Período (meses) usados na geração de uma base de suporte (img. mediana)
                               sensorName,                #Nome do sensor (útil a criação das máscaras de nuvem)
                               cloudCoverMaxPercentage,   #Percentual de cobertura de nuvem
                               UsefulnessThres):          #Limiar de aproveitamento da imagem (cobertura de nuvem+background)

    radiusMask = 2 

    #Obter informações do sensor
    info = specIndex.sensor_info(sensorName['multispecHR'])
    BoI = info['BoI']

    repack = []
    listOfInstants = []
    dim = collection.size().getInfo()
    imageList = collection.toList(collection.size())
    for i in range(dim):
        imageInfo = ee.Image(imageList.get(i)).getInfo()
        image = ee.Image(imageList.get(i))

        image = ee.Image(-1).blend(image) #Desnecessário?
        #Definir a máscara desta imagem (melhorar a parametrização)
        cloudMask = build_cloud_mask_V2(image, sensorName, 'multispecHR',geometry,scale)

        #Definição da máscara de nuvem+background (i.e., dados faltante)
        cloudBackgroundMask_out = ee.Image(0).blend(cloudMask).rename('cloudAndBackgroundMask')

        #Com base no comportamento de "cloudBackgroundMask_out", é possível ter uma noção do aproveitamento da imagem!
        #Obtém e transforma a máscara e uma lista 0-1 que indica a pertinencia ou não em área útil
        valuesMask = np.array( ee.List(cloudBackgroundMask_out \
                                 .reduceRegion(reducer=ee.Reducer.toList(),geometry=geometry,scale=scale,bestEffort=True)) \
                                 .getInfo()['cloudAndBackgroundMask'] )
        #Realiza o cálculo do percentual de elementos "1" (i.e., em área útil) em relação ao total
        Usefulness = np.count_nonzero(valuesMask == 1) / valuesMask.shape[0]
        
        if Usefulness >= UsefulnessThres:

            #------------------------------------------------
            #Definir uma mediana em relação ao instante de "imageList.get(i))
            #imagingTime = imageInfo['properties']['SENSING_TIME'].split('T')[0]
            imagingTime =  (imageInfo['properties']['DATATAKE_IDENTIFIER'].split('_')[1]).split('T')[0]

            #imagingTime = imageInfo['properties']['DATE_ACQUIRED']
            #endYMD = datetime.strptime(imagingTime, "%Y-%m-%d")     #mais recente (fim)      
            endYMD = datetime.strptime(imagingTime, "%Y%m%d")     #mais recente (fim)      

            ref = endYMD - relativedelta.relativedelta(months=pastMonthsMedian)  #3 é um teste...
            backTime = str(ref.year)+'-'+str(ref.month)+'-'+str(ref.day)
            endTime = str(endYMD.year)+'-'+str(endYMD.month)+'-'+str(endYMD.day)
            #Definir a busca de modo mais conveniente aqui (melhorar a parametrização)
            earlyCollection = ee.ImageCollection(info['completeSensorName']) \
                                .filterBounds(geometry) \
                                .filterDate(backTime,endTime) \
                                .filterMetadata('CLOUDY_PIXEL_PERCENTAGE','less_than',cloudCoverMaxPercentage)

            #earlyCollection = collection_chrono_sort(earlyCollection) #<<<Desnecessário reordenar...

            medianImage = earlyCollection.median()
            #------------------------------------------------

            #Aplicar um filtro de dilatação (circular de raio 1)
            cloudBackgroundMask_out = cloudBackgroundMask_out.focal_max(kernel=ee.Kernel.circle(radius=radiusMask), iterations=1)

            #Adicao da mascara de dados faltantes como banda da imagem (para fins de referencia futura)
            image = image.addBands([cloudBackgroundMask_out])

            #Aplicação das máscaras In/Out sobre as BoI's
            newBandNames = []
            for band in BoI:

                #Imagem/banda original -- aplicação da máscara para remoção dos dados ausentes/contaminados
                __bandMasked_in  = image.select(band).updateMask(cloudBackgroundMask_out).rename(band+'_mask'+'In')
            
                #Low-level median -- aplicação da máscara reversa para aproveitamento dos dados
                __bandMedianMasked_out = medianImage.select(band).updateMask(cloudBackgroundMask_out.Not()).rename(band+'median_mask'+'Out')

                #Preenchimento de dados faltantes
                __blendBand = __bandMedianMasked_out.blend(__bandMasked_in).rename(band+'_corr')

                image = image.addBands([__blendBand])

            #Renomear as bandas corrigidas e originais
            bandNames = image.bandNames().getInfo()
            #for b in bandNames:
            for b in BoI:
                posOrig = bandNames.index(b)
                posCorr = bandNames.index(b+'_corr')
                bandNames[posCorr] = bandNames[posOrig]
                bandNames[posOrig] = bandNames[posOrig]+'_original' 
            image = image.select(image.bandNames().getInfo()).rename(bandNames)

            #Armazenando a imagem para instanciar a coleção em seguida
            repack.append(image)

            #Armazenar o instante da imagem extraída da coleção
            listOfInstants.append( str(datetime.strptime(imagingTime, "%Y%m%d")).split(' ')[0] )

            #Aplicar um coletor de lixo aqui?
            gc.collect()

    #Reempacotamento da coleção, após ajuste de nuvens/background
    repackCollection = ee.ImageCollection.fromImages( repack )
    
    return repackCollection, listOfInstants

#======================================
#Versão preliminar que só funciona para Landsat
def collection_chrono_sort(collection):#,sensorName):

    dim = collection.size().getInfo()
    imageList = collection.toList(collection.size())
    vecTime = []
    for i in range(dim):
        imageInfo = ee.Image(imageList.get(i)).getInfo()
        temp =  datetime.strptime(imageInfo['properties']['SENSING_TIME'].split('T')[0], "%Y-%m-%d")
        vecTime.append( temp )
    #sortVec = np.sort(vecTime)
    argsortVec = np.argsort(vecTime)
    
    #Reempacotando a coleção ordenada...
    repack = []
    for index in argsortVec:
        image = ee.Image(imageList.get(int(index)))
        repack.append( image )
    sortCollection = ee.ImageCollection.fromImages( repack )

    return sortCollection, vecTime

#======================================
def collection_chrono_sort_sentinel(collection):#,sensorName):

    dim = collection.size().getInfo()
    imageList = collection.toList(collection.size())
    vecTime = []
    for i in range(dim):
        imageInfo = ee.Image(imageList.get(i)).getInfo()
        temp = datetime.strptime( (imageInfo['properties']['DATATAKE_IDENTIFIER'].split('_')[1]).split('T')[0] , '%Y%m%d')
        vecTime.append( temp )
    #sortVec = np.sort(vecTime)
    argsortVec = np.argsort(vecTime)
    
    #Reempacotando a coleção ordenada...
    repack = []
    for index in argsortVec:
        image = ee.Image(imageList.get(int(index)))
        repack.append( image )
    sortCollection = ee.ImageCollection.fromImages( repack )

    return sortCollection

#======================================
#Versão preliminar que só funciona para Landsat
def collection_chrono_sort__25set21(collection):#,sensorName):

    dim = collection.size().getInfo()
    imageList = collection.toList(collection.size())
    vecTime = []
    for i in range(dim):
        imageInfo = ee.Image(imageList.get(i)).getInfo()
        temp =  datetime.strptime(imageInfo['properties']['SENSING_TIME'].split('T')[0], "%Y-%m-%d")
        vecTime.append( temp )
    #sortVec = np.sort(vecTime)
    argsortVec = np.argsort(vecTime)
    
    #Reempacotando a coleção ordenada...
    repack = []
    chronoPack = []
    for index in argsortVec:
        image = ee.Image(imageList.get(int(index)))
        repack.append( image )
        chronoPack.append(vecTime[index])
    sortCollection = ee.ImageCollection.fromImages( repack )

    return sortCollection, chronoPack

#--------------------------------------------
#Função para análise e adequação da coleção
def cloud_shadow_fill(image,                #Coleção sob análise
                      startDate,
                      endDate,
                      geometry,                  #Geometria da coleção
                      scale,                     #Escala dos dados
                      pastMonthsMedian,          #Período (meses) usados na geração de uma base de suporte (img. mediana)
                      sensorName,                #Nome do sensor (útil a criação das máscaras de nuvem)
                      cloudCoverMaxPercentage):  #Percentual de cobertura de nuvem                      

    radiusMask = 2 #estava 1... deve melhorar com 2... (?)

    #Obter informações do sensor
    info = specIndex.sensor_info(sensorName['multispec'])
    BoI = info['BoI']


    imageInfo = image.getInfo()
    image = ee.Image(-1).blend(image) 
        
    #Definir a máscara desta imagem (melhorar a parametrização)
    cloudMask = build_cloud_mask_V2(image, sensorName, 'multispec',geometry,scale)

    #Definição da máscara de nuvem+background (i.e., dados faltante)
    cloudBackgroundMask_out = ee.Image(0).blend(cloudMask).rename('cloudAndBackgroundMask')

    #Com base no comportamento de "cloudBackgroundMask_out", é possível ter uma noção do aproveitamento da imagem!
    #Obtém e transforma a máscara e uma lista 0-1 que indica a pertinencia ou não em área útil
    valuesMask = np.array( ee.List(cloudBackgroundMask_out \
                                 .reduceRegion(reducer=ee.Reducer.toList(),geometry=geometry,scale=scale,bestEffort=True)) \
                                 .getInfo()['cloudAndBackgroundMask'] )

    #------------------------------------------------
    #Definir uma mediana em relação ao instante de "imageList.get(i))
    imagingTime = imageInfo['properties']['SENSING_TIME'].split('T')[0]
    
    endYMD = datetime.strptime(imagingTime, "%Y-%m-%d")     #mais recente (fim)      

    ref = endYMD - relativedelta.relativedelta(months=pastMonthsMedian)  #3 é um teste...
    backTime = str(ref.year)+'-'+str(ref.month)+'-'+str(ref.day)
    endTime = str(endYMD.year)+'-'+str(endYMD.month)+'-'+str(endYMD.day)
            
    # #Definir a busca de modo mais conveniente aqui (melhorar a parametrização)
    # earlyCollection = ee.ImageCollection(info['completeSensorName']) \
    #                             .filterBounds(geometry) \
    #                             .filterDate(backTime,endTime) \
    #                             .filterMetadata('CLOUDY_PIXEL_PERCENTAGE','less_than',cloudCoverMaxPercentage)    
    earlyCollection = ee.ImageCollection(info['completeSensorName']) \
                       .filterBounds(geometry) \
                       .filterDate(backTime,endTime) \
                       .filterMetadata('CLOUD_COVER','less_than', cloudCoverMaxPercentage)

    medianImage = earlyCollection.median()
    #------------------------------------------------

    #Aplicar um filtro de dilatação (circular de raio 1)
    cloudBackgroundMask_out = cloudBackgroundMask_out.focal_max(kernel=ee.Kernel.circle(radius=radiusMask), iterations=1)

    #Adicao da mascara de dados faltantes como banda da imagem (para fins de referencia futura)
    image = image.addBands([cloudBackgroundMask_out])

    #Aplicação das máscaras In/Out sobre as BoI's
    #newBandNames = []
    for band in BoI:

        #Imagem/banda original -- aplicação da máscara para remoção dos dados ausentes/contaminados
        __bandMasked_in  = image.select(band).updateMask(cloudBackgroundMask_out).rename(band+'_mask'+'In')
            
        #Low-level median -- aplicação da máscara reversa para aproveitamento dos dados
        __bandMedianMasked_out = medianImage.select(band).updateMask(cloudBackgroundMask_out.Not()).rename(band+'median_mask'+'Out')

        #Preenchimento de dados faltantes
        __blendBand = __bandMedianMasked_out.blend(__bandMasked_in).rename(band+'_corr')

        image = image.addBands([__blendBand])

    #Renomear as bandas corrigidas e originais
    #bandNames = image.bandNames().getInfo()
    #for b in bandNames:
    for b in BoI:
        bandNames = image.bandNames().getInfo()
        posOrig = bandNames.index(b)
        posCorr = bandNames.index(b+'_corr')
        bandNames[posCorr] = bandNames[posOrig]
        bandNames[posOrig] = bandNames[posOrig]+'_original' 
        image = image.select(image.bandNames().getInfo()).rename(bandNames)

    #Aplicar um coletor de lixo aqui?
    gc.collect()

    return image
















