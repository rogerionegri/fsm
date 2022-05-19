#--------------------------------------------
#IMPORTS ------------------------------------
import numpy as np
import pandas as pd
from datetime import datetime
from dateutil.relativedelta import relativedelta
import gc #garbage collector
from sklearn import preprocessing

#--------------------------------------------
#GEE ----------------------------------------
import ee
ee.Initialize() 

#--------------------------------------------
#PERSONAL MODULES ---------------------------
import specIndex
import fireAnalysis
import maskingFunctions

#--------------------------------------------
def get_geometry_from_corners(coordenadas):
    
    x1,y1,x2,y2 = coordenadas.split(",")
    geometria = ee.Geometry.Polygon( \
                                    [[[float(x1),float(y2)], \
                                      [float(x2),float(y2)], \
                                      [float(x2),float(y1)], \
                                      [float(x1),float(y1)], \
                                      [float(x1),float(y2)]]])

    return geometria

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
#Adotar um esquema recursivo sobre "np" (que retorna a data cada vez mais, até que encontre alguma imagem)
#Contempla apenas indices espectrais 
def expand_period_collection(refData,T,geometry,scale,spamMedianSupport,UsefulnessThres,sensorName,cloudCoverMaxPercentage,coordSystem):
    
    print("Busca de nivel "+str(T))

    info = specIndex.sensor_info(sensorName['multispec'])

    ymd = datetime.strptime(refData, "%Y-%m-%d")    
    startYDM = ymd - relativedelta(months=T)
    startData = str(startYDM.year)+'-'+str(startYDM.month)+'-'+str(startYDM.day)
    endData = refData
    
    pre_collection = ee.ImageCollection(info['completeSensorName']) \
                       .filterBounds(geometry) \
                       .filterDate(startData,endData) \
                       .filterMetadata('CLOUD_COVER','less_than', cloudCoverMaxPercentage)

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
#Contempla apenas indices espectrais 
def get_df_inst(refData,
                nPast,
                geometry,
                scale,
                sensorName,
                cloudCoverMaxPercentage,
                coordSystem,
                getFocos):

    defaultDummy = 0.0 #-9999.0

    #Obtenha informações básicas do sensor adotado
    infoMSpec = specIndex.sensor_info(sensorName['multispec'])
    infoPrec = specIndex.sensor_info(sensorName['precipitation'])
    infoTemp = specIndex.sensor_info(sensorName['temperature'])
    infoHuEv = specIndex.sensor_info(sensorName['humidity'])     
    infoFire = specIndex.sensor_info(sensorName['fire'])
    
    #Define a data type useful to manipulate time variable
    ymd = datetime.strptime(refData, "%Y-%m-%d")

    d = {} #Definindo uma estrutura de dicionário para construir um dataFrame no final do processo
    
    #Faça um loop no intervalo de tempo "nPast" (em meses)
    for t in range(1,nPast+1):
        
        #Data final para a coleção do mês 
        endYDM = ymd - relativedelta(months=t-1)
        endData = str(endYDM.year)+'-'+str(endYDM.month)+'-'+str(endYDM.day)

        #Util no caso dos dados de precipitacao e temperatura
        startYDM = ymd - relativedelta(months=t)
        startData = str(startYDM.year)+'-'+str(startYDM.month)+'-'+str(startYDM.day)

        #Dados multiespectrais-------------------------------
        monthCollection = expand_period_collection(endData,1,geometry,scale,sensorName,cloudCoverMaxPercentage,coordSystem)

        #Cálculo de índices espectrais
        monthCollection = specIndex.add_ndvi_collection(monthCollection,sensorName['multispec'])
        monthCollection = specIndex.add_ndwi_collection(monthCollection,sensorName['multispec'])
        monthCollection = specIndex.add_savi_collection(monthCollection,sensorName['multispec'])
        monthCollection = specIndex.add_nbr_collection(monthCollection,sensorName['multispec'])

        #Cálculo da "mediana espectral" no mês
        medianMonth = monthCollection.median()

        #Extração dos dados sobre índices espectrais
        #Utilizando lonMS/latMS em distinção aos demais para fins de teste (e inclusão no dict ao final do processo)
        temp = ee.Image(defaultDummy).blend(ee.Image( medianMonth ))
        latMS, lonMS, indNDVI = ext_lat_lon_pixel(temp, geometry, ['ndvi'], scale, coordSystem)
        latMS, lonMS, indNDWI = ext_lat_lon_pixel(temp, geometry, ['ndwi'], scale, coordSystem)
        latMS, lonMS, indSAVI = ext_lat_lon_pixel(temp, geometry, ['savi'], scale, coordSystem)  
        latMS, lonMS, indNBR  = ext_lat_lon_pixel(temp, geometry, ['nbr'], scale, coordSystem)   

        #Dados de precipitação-------------------------------
        monthPrecipitation = ee.ImageCollection(infoPrec['completeSensorName']) \
                               .filterBounds(geometry) \
                               .filterDate(startData,endData)

        #Cálculo da precipitação acumulada no mês
        cumMonthPrec = monthPrecipitation.sum().rename('precipitation') #<<<funciona pq tem uma banda só!!
        
        #Fazer a adição da precipitação do mês
        temp = temp.addBands( ee.Image(cumMonthPrec) )

        #Extração dos dados sobre preciptação acumulada
        lat, lon, indPrec = ext_lat_lon_pixel(temp, geometry, ['precipitation'], scale, coordSystem)

        #Dados de temperatura----------------------------------
        monthTemperature = ee.ImageCollection(infoTemp['completeSensorName']) \
                             .filterBounds(geometry) \
                             .filterDate(startData,endData)

        #Cálculo da temperatura média no mês
        #avgMonthTemp = monthTemperature.mean().rename('temperature')
        #"Select" apenas na banda de interesse... ('temperature_2m')
        avgMonthTemp = monthTemperature.mean().select(infoTemp['bands']['temperature_2m']).rename('temperature')
        
        #Fazer a adição da temperatura do mês
        temp = temp.addBands( ee.Image(avgMonthTemp) )

        #Extração dos dados sobre preciptação acumulada
        lat, lon, indTemp = ext_lat_lon_pixel(temp, geometry, ['temperature'], scale, coordSystem)

        #Dados de evapotranspiração e umidade----------------------
        monthHumidtyEvapo = ee.ImageCollection(infoHuEv['completeSensorName']) \
                              .filterBounds(geometry) \
                              .filterDate(startData,endData)

        #Cálculo da evapotranspiração e humidade média no mês
        avgMonthHuEv = monthHumidtyEvapo.mean().select( \
                                                         [ infoHuEv['bands']['specificHumidity'] ,   \
                                                           infoHuEv['bands']['Evapotranspiration'] ] \
                                                      ).rename(['humidity','evapotransp'])

        #Fazer a adição da umidade e evapotranspiração do mês
        temp = temp.addBands( ee.Image(avgMonthHuEv) )

        #Extração dos dados sobre umidade e evapotranspiração média do mês
        lat, lon, indHuEv = ext_lat_lon_pixel(temp, geometry, ['humidity','evapotransp'], scale, coordSystem)

        #>>> SE A KEYWORD ESTIVER ATIVADA!
        if getFocos:
            #Dados do produto FIRMS---------------------------------
            monthFire = ee.ImageCollection(infoFire['completeSensorName']) \
                          .filterBounds(geometry) \
                          .filterDate(startData,endData)

            #lon, lat, indNDVI = ext_lat_lon_pixel(temp, geometry, ['ndvi'], scale, coordSystem)
            monthFire = specIndex.add_firespots_collection(monthFire,sensorName['fire'])

            countMonthFires = monthFire.sum().select('occurrence').rename('contOccurrence')
            #infoTemp['bands']['temperature_2m']

            temp = temp.addBands( ee.Image(countMonthFires) )
            #Extração dos dados sobre os focos de incêndio
            lat, lon, indFocos = ext_lat_lon_pixel(temp, geometry, ['contOccurrence'], scale, coordSystem)
       
        #-----------------------------------------------
        d['NDVI'+str(t)] = indNDVI[0]    #<<< [0] for compatibility (this is a numpyArray)
        d['NDWI'+str(t)] = indNDWI[0]
        
        d['SAVI'+str(t)] = indSAVI[0]   #Add 19mar21
        d['NBR'+str(t)] = indNBR[0]     #Add 19mar21
        
        d['Prec'+str(t)] = indPrec[0]
        d['Temp'+str(t)] = indTemp[0]    #<<< estou achando estranho os valores, que deveriam estar em kelvin
        
        d['Humi'+str(t)] = indHuEv[0] * (10**2)   #Add 19mar21 (uso de fator mutiplicativo)
        d['Evap'+str(t)] = indHuEv[1] * (10**5)   #Add 19mar21

        #Armazenar a informação sobre a ocorrência de focos
        if getFocos: d['Foco'+str(t)] = indFocos[0]

        print('>>>'+str(t))

        gc.collect() #Coletor de lixo :: Precisa desalocar a var. "temp"?

    #Inclusão da lat/lon no dict
    d['Latitude'] = latMS
    d['Longitude'] = lonMS

    #Construct and fill in the resulting dataFrame
    tab = pd.DataFrame()
    tab = tab.from_dict(d)   #<<<Algumas variáveis (prec e temp) não preenchem toda a geometria... isso pode causar probelma aqui!


    if getFocos:
        #tabRev = fireAnalysis.fire_spots_df_v2(tab,refData)#,gmmComponents)
        #tabRev = fireAnalysis.fire_spots_df_v3(tab,refData)
        tabRev = fireAnalysis.fire_spots_df_v4(tab,refData)
        return tabRev
    else:
        return tab

#-----------------------------------------------
#Organiza dos dados entre treino e predição, além de separar os atributos de interesse e normaliza-los
def train_predict_data_frames(__trainBase,__predictBase,listVars):

    #Cópias de segurança (verificar se há necessidade de fato...)
    _trainBase = __trainBase.copy()
    _predictBase = __predictBase.copy()

    #Inicialização de dicionarios para definir as bases normalizadas
    dTrain = {}
    dPredict = {}

    #Lista de atributos
    listAtts = []

    #Instanciação do método de normalizaçao
    normalizador = preprocessing.MinMaxScaler()
    for item in listVars:
    
        coi = [col for col in _trainBase.columns if item in col]  #...column of interest
        listAtts = listAtts + coi

        for var in coi: #...já é o nome das colunas na forma correta

            xT = _trainBase.loc[:,var].values.astype(float)
            xP = _predictBase.loc[:,var].values.astype(float)  #os nomes entre tabelas devem estar iguais, pois foram gerados pela mesma função

            #Definição da função normalizadora
            normalizador.fit( (xT.T).reshape(-1, 1) )

            #Normalização dos dados de treinamento
            xT_scaled = normalizador.transform( xT.T.reshape(-1,1) )
            #Normalização dos dados para predição
            xP_scaled = normalizador.transform( xP.T.reshape(-1,1) )

            dTrain[var] = xT_scaled[:,0]
            dPredict[var] = xP_scaled[:,0]

    #Montagem das bases
    dTrain['Longitude'] = _trainBase.Longitude
    dTrain['Latitude'] = _trainBase.Latitude
    dTrain['prob'] = _trainBase.prob
    dTrain['logProb'] = _trainBase.logProb

    dPredict['Longitude'] = _predictBase.Longitude
    dPredict['Latitude'] = _predictBase.Latitude

    trainBase = pd.DataFrame(dTrain)
    predictBase = pd.DataFrame(dPredict)

    print("Pit stop...")

    #Dataframes...
    #...padrões
    X = trainBase.loc[:,listAtts].values.astype(float)
    #...rótulos
    Y = trainBase.loc[:,'prob'].values.astype(float)
    #...padrões não rotulados para predição
    P = predictBase.loc[:,listAtts].values.astype(float)

    return X, Y, P, trainBase, predictBase

#-----------------------------------------------
#Organiza dos dados entre treino e predição, além de separar os atributos de interesse e normaliza-los
def train_predict_data_frames_v2(__trainBase,__predictBase,listVars):

    #Cópias de segurança (verificar se há necessidade de fato...)
    _trainBase = __trainBase.copy()
    _predictBase = __predictBase.copy()

    #Inicialização de dicionarios para definir as bases normalizadas
    dTrain = {}
    dPredict = {}

    #Lista de atributos
    listAtts = []

    #Instanciação do método de normalizaçao
    normalizador = preprocessing.MinMaxScaler()
    for item in listVars:
    
        coi = [col for col in _trainBase.columns if item in col]  #...column of interest
        listAtts = listAtts + coi

        for var in coi: #...já é o nome das colunas na forma correta

            xT = _trainBase.loc[:,var].values.astype(float)
            xP = _predictBase.loc[:,var].values.astype(float)  #os nomes entre tabelas devem estar iguais, pois foram gerados pela mesma função

            #Definição da função normalizadora (sobre os dados de PREDICAO -- mais genérico/geral)
            #normalizador.fit( (xT.T).reshape(-1, 1) )
            normalizador.fit( (xP.T).reshape(-1, 1) )

            #Normalização dos dados de treinamento
            xT_scaled = normalizador.transform( xT.T.reshape(-1,1) )
            #Normalização dos dados para predição
            xP_scaled = normalizador.transform( xP.T.reshape(-1,1) )

            dTrain[var] = xT_scaled[:,0]
            dPredict[var] = xP_scaled[:,0]

    #Montagem das bases
    dTrain['Longitude'] = _trainBase.Longitude
    dTrain['Latitude'] = _trainBase.Latitude
    #dTrain['indicador'] = _trainBase.indicador  #Talvez seja até desnecessário...
    #dTrain['logProb'] = _trainBase.logProb

    dPredict['Longitude'] = _predictBase.Longitude
    dPredict['Latitude'] = _predictBase.Latitude

    trainBase = pd.DataFrame(dTrain)
    predictBase = pd.DataFrame(dPredict)

    print("Pit stop...")


    #Dataframes...
    #...padrões
    X = trainBase.loc[:,listAtts].values.astype(float)
    #...rótulos
    #Y = trainBase.loc[:,'prob'].values.astype(float)
    #...padrões não rotulados para predição
    P = predictBase.loc[:,listAtts].values.astype(float)

    #return X, Y, P, trainBase, predictBase
    return X, P, trainBase, predictBase

#-----------------------------------------------
#Organiza dos dados entre treino e predição, além de separar os atributos de interesse e normaliza-los
def train_predict_data_frames_v3(__trainBase,__predictBase,listVars):

    #Cópias de segurança 
    _trainBase = __trainBase.copy()
    _predictBase = __predictBase.copy()

    #Inicialização de dicionarios para definir as bases normalizadas
    dTrain = {}
    dPredict = {}

    #Lista de atributos
    listAtts = []

    #Instanciação do método de normalizaçao
    normalizador = preprocessing.MinMaxScaler()
    for item in listVars:
    
        coi = [col for col in _trainBase.columns if item in col]  #...column of interest
        listAtts = listAtts + coi

        for var in coi: #...já é o nome das colunas na forma correta

            xT = _trainBase.loc[:,var].values.astype(float)
            xP = _predictBase.loc[:,var].values.astype(float)  #os nomes entre tabelas devem estar iguais, pois foram gerados pela mesma função

            temp = np.concatenate([xP,xT]) #xT.append(xP)

            #Definição da função normalizadora 
            normalizador.fit( (temp.T).reshape(-1, 1) )

            #Normalização dos dados de treinamento
            xT_scaled = normalizador.transform( xT.T.reshape(-1,1) )

            #Normalização dos dados para predição
            xP_scaled = normalizador.transform( xP.T.reshape(-1,1) )

            dTrain[var] = xT_scaled[:,0]
            dPredict[var] = xP_scaled[:,0]

    #Montagem das bases
    dTrain['Longitude'] = _trainBase.Longitude
    dTrain['Latitude'] = _trainBase.Latitude

    dPredict['Longitude'] = _predictBase.Longitude
    dPredict['Latitude'] = _predictBase.Latitude

    trainBase = pd.DataFrame(dTrain)
    predictBase = pd.DataFrame(dPredict)

    print("Pit stop...")

    #Dataframes...
    #...padrões
    X = trainBase.loc[:,listAtts].values.astype(float)
    #...rótulos
    P = predictBase.loc[:,listAtts].values.astype(float)

    #return X, Y, P, trainBase, predictBase
    return X, P, trainBase, predictBase

#-----------------------------------------------
#Organiza dos dados entre treino e predição, além de separar os atributos de interesse e normaliza-los
def train_predict_data_frames_v4(__trainBase,__predictBase,listVars):

    #Cópias de segurança (verificar se há necessidade de fato...)
    _trainBase = __trainBase.copy()
    _predictBase = __predictBase.copy()

    #Inicialização de dicionarios para definir as bases normalizadas
    dTrain = {}
    dPredict = {}

    #Lista de atributos
    listAtts = []

    #Instanciação do método de normalizaçao
    normalizador = preprocessing.MinMaxScaler()
    for item in listVars:
    
        coi = [col for col in _trainBase.columns if item in col]  #...column of interest
        listAtts = listAtts + coi

        xT = _trainBase.loc[:,coi].values.astype(float).max(axis=1)
        xP = _predictBase.loc[:,coi].values.astype(float).max(axis=1)
        
        temp = np.concatenate([xP,xT])

        normalizador.fit( (temp.T).reshape(-1, 1) )

        #Normalização dos dados de treinamento
        xT_scaled = normalizador.transform( xT.T.reshape(-1,1) )
        #Normalização dos dados para predição
        xP_scaled = normalizador.transform( xP.T.reshape(-1,1) )

        dTrain[item] = xT_scaled[:,0]
        dPredict[item] = xP_scaled[:,0]

    #Montagem das bases
    dTrain['Longitude'] = _trainBase.Longitude
    dTrain['Latitude'] = _trainBase.Latitude

    dPredict['Longitude'] = _predictBase.Longitude
    dPredict['Latitude'] = _predictBase.Latitude

    trainBase = pd.DataFrame(dTrain)
    predictBase = pd.DataFrame(dPredict)

    print("Pit stop...")

    #Dataframes...
    #...padrões
    X = trainBase.loc[:,listVars].values.astype(float)
    #...padrões não rotulados para predição
    P = predictBase.loc[:,listVars].values.astype(float)

    return X, P, trainBase, predictBase