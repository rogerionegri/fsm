#--------------------------------------------
#GEE ----------------------------------------
import ee
ee.Initialize() 


def get_modis_fire_image(startTime,endTime,geometry,scale,coordSys):

    #Calculando índice 'focos' para ocorrência de fogo
    def spots_modis(image):
        focos = image.expression('int((BurnDate > 0) && (Uncertainty < 10))', {'BurnDate': image.select('BurnDate'), 'Uncertainty': image.select('Uncertainty')}).rename('focos')
        return image.addBands([focos])


    #Calculando índice 'focos' para ocorrência de fogo
    def spots_modis_v2(image):
        valBurnData, valUncertity = 0, 50
        focos = image.expression('int((BurnDate > valBD) && (Uncertainty < valU))', {'valBD': valBurnData, 'valU': valUncertity, 'BurnDate': image.select('BurnDate'), 'Uncertainty': image.select('Uncertainty')}).rename('focos')
        return image.addBands([focos])


    #=======================================
    #Definições internas...
    #Banda de interesse: focos 
    bandInteresse = ['focos'] 
    scale = 100 # --> a escala do MODIS é de 500m
    defaultDummy = 0 #-9990.0
    restricted = 1
    #=======================================

    #Consultando coleção: MCD64A1.006 MODIS Burned Area Monthly Global 500m 
    collection = ee.ImageCollection('MODIS/006/MCD64A1') \
                   .filterBounds(geometry) \
                   .filterDate(startTime,endTime) \
                   .sort('system:time_start', True)
                   
    print("Total de imagens encontradas: "+str(collection.size().getInfo()))
    if restricted: 
        collection = collection.map(spots_modis) 
    else:
        collection = collection.map(spots_modis_v2)

    #=======================================
    #statImg = collection.sum()
    statImg = collection.sum().unmask(defaultDummy)
    refImage = ee.Image(defaultDummy).blend(statImg).clipToBoundsAndScale(geometry=geometry,scale=scale)
    
    return refImage
