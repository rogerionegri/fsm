#--------------------------------------------
#IMPORTS ------------------------------------
import pandas as pd
import numpy as np

#--------------------------------------------
#GEE ----------------------------------------
import ee
ee.Initialize() 

#--------------------------------------------
#PERSONAL MODULES ---------------------------
import dataManage
import regressionModels
import saveModule
import configModule
import MFM_DataBase
import UFD_DataBase  
import moduleUtils   
import moduleModis

#---------------------------------------------
#PARAMETERS ----------------------------------
#---------------------------------------------

#Basic configurations-------------------------
coordSystem = "epsg:4326"
scale = 100  
cloudCoverMaxPercentage = 50.0  
spamMedianSupport = 3  
UsefulnessThres = 0.75 
statTypeRef = 'ref'  
minFireOccur = 1  

#---------------------------------------------------
#Additional parameters (not change)
listAtts = ['NDVI','NDWI','NBR']
alpha = 0.1     
beta  = 0.0
early_fire = 3   
                 
                 
#Anomaly detection model: 'IsolationForest' or 'OCSVM'
nJobs = 4                       
modelName = 'IsolationForest'   

#Study area geometry --------------------------------
coordinates = "-58.01999969470172,-16.57724279945931,-56.99003143298297,-16.026287800090003" #Example: Horizonte d'Oeste area 

#Pre-fire period
spamPeriod_init = "2015-01-01" #Example of periods
spamPeriod_end  = "2017-01-01" #Please, change according your convenience

#Modeling period
analysisPeriod_init = "2017-06-01"
analysisPeriod_end  = "2018-03-31"

#Assessment period
assessmentPeriod_init = "2018-07-01"
assessmentPeriod_end = "2018-08-31"

#Mininum instant lapse
lapse = 15 #in days

#Output path
path_out = './'

#Option vector (Yes/No) --> save preliminar results
#[imgRef,reports,modisSuport,modisRef,baseImage,outputInstantMaps,outputMaps]  
optionSaves = [0,1,0,1,1,1,1]   
#----------------------------------------------

#Init of processings---------------------------
sensorName = configModule.get_sensor_name()

#Defining a "geometry" object regarding the study area/region of interest
geometry = dataManage.get_geometry_from_corners(coordinates)


#------------------------------------------
#Compute the pre-fire image----------------------------------------
boolFlag, baselinePrefireCollection, LoI_Prefire = UFD_DataBase.baseline_prefire_norec(spamPeriod_init,\
                                                          spamPeriod_end,\
                                                          geometry,\
                                                          scale,\
                                                          spamMedianSupport,\
                                                          UsefulnessThres,\
                                                          sensorName,\
                                                          cloudCoverMaxPercentage,\
                                                          coordSystem)

prefireMedian, baselinePrefireCollection = UFD_DataBase.include_attributes(baselinePrefireCollection, \
                                                                           sensorName, \
                                                                           listAtts,\
                                                                           'baseline', \
                                                                           'median')

#Save the reference/trend image?
if optionSaves[0]:
    moduleUtils.save_image(prefireMedian, path_out+'refImage.tif', geometry, scale, coordSystem)
    moduleUtils.add_instant_info(path_out+'refImage',LoI_Prefire)


#------------------------------------------
#Compute the "post-fire" images (ref. modeling period)
flag, postfireSerie, LoI = UFD_DataBase.timeseries_postfire_simple(analysisPeriod_init,
                                                 analysisPeriod_end,
                                                 listAtts,
                                                 geometry,
                                                 scale,
                                                 spamMedianSupport,
                                                 UsefulnessThres,
                                                 sensorName,
                                                 cloudCoverMaxPercentage,
                                                 coordSystem)

#Save the "post-fire"/modeling images?
if optionSaves[1]:
    moduleUtils.add_report_info(path_out+'report.txt','Modeling period image',LoI)

#------------------------------------------
baseModis = moduleModis.get_modis_fire_image(analysisPeriod_init,analysisPeriod_end,geometry,scale,coordSystem)
#baseModisRef = moduleModis.get_modis_fire_image(modisAnalysisPeriod_init,modisAnalysisPeriod_end,geometry,scale,coordSystem)

#Save the Modis image series used as anciliary data to identify fire events?
if optionSaves[2]:
    moduleUtils.save_image(baseModis, path_out+'modisSupport.tif', geometry, scale, coordSystem)

#------------------------------------------
#Compute the deltaNBR
baseImage = MFM_DataBase.build_base_image__varModis(prefireMedian,
                                          postfireSerie,
                                          listAtts,
                                          'baseline',
                                          'mean',
                                          alpha,
                                          beta,
                                          early_fire,
                                          geometry,
                                          scale,
                                          coordSystem,
                                          baseModis)

    
#Save the pos-fire/deltaNBR images
if optionSaves[4]:
    moduleUtils.save_image(baseImage, path_out+'baseImage.tif', geometry, scale, coordSystem)
    
#------------------------------------------
#Defining a dataFrame with extracted data
baseDataframe = MFM_DataBase.build_base_dataframe(baseImage,
                                     listAtts,
                                     'ref', 
                                     geometry,
                                     scale,
                                     coordSystem)

#==========================================================================
#Susceptibility mapping (assessment period)---------------
spamMedianSupport__ = 2 #non-adjustable parameter
                        
assessPeriod = MFM_DataBase.early_assessment_period(assessmentPeriod_init,\
                                                    assessmentPeriod_end,\
                                                    lapse,\
                                                    listAtts,\
                                                    geometry,\
                                                    scale,\
                                                    spamMedianSupport__,\
                                                    UsefulnessThres,\
                                                    sensorName,\
                                                    cloudCoverMaxPercentage,\
                                                    coordSystem)

assessmentImage, dictImageInstants, cloudMasksInstants = MFM_DataBase.build_assessment_early_instants_image(
                                                      assessPeriod,
                                                      listAtts)#,#geometry,#sensorName,#coordSystem)

#Apply the anomaly detection model on each instante of assessmentPeriod_init~end
detectorModel = regressionModels.get_detector(modelName,nJobs)
data = ( baseDataframe.loc[ baseDataframe.sum_thres_deltaNBR > minFireOccur , [i+'_'+statTypeRef  for i in listAtts] ] ).to_numpy()
detectorModel.fit(data)    

#Defining a dataFrame with assessment data
dictInstMaps = {}
for instant in dictImageInstants:
    
    instAssessDataframe = MFM_DataBase.build_assessment_instant_dataframe__withClouds(assessmentImage,
                                                          #dictImageInstants,
                                                          cloudMasksInstants,
                                                          int(instant),
                                                          listAtts,
                                                          geometry,
                                                          scale,
                                                          coordSystem)#,sensorName)

    
    instData = ( instAssessDataframe.loc[ : , listAtts ] ).to_numpy()
    print('Predict in... '+str(instant)+'/'+str(len(dictImageInstants)))
    instMap = detectorModel.predict(instData)
    print('...predict out '+str(instant))
    
    instAssessDataframe['res'] = instMap
    instAssessDataframe.loc[ (instAssessDataframe.loc[:,'cloud'] != 1) , 'res'] = 0 #locations with cloud in some instant will have mean value dfferent of 1

    #Feeding a dictionary...
    dictInstMaps['map_'+str(instant)] = instAssessDataframe.res.to_numpy()
    dictInstMaps['cloud_'+str(instant)] = ((instAssessDataframe['cloud']).to_numpy() == 1) #Test where does not occur cloud/shadow
    dictInstMaps['anomaly_'+str(instant)] = (dictInstMaps['map_'+str(instant)] < 0)

    #Save the images in the analysis period
    if optionSaves[5]:
        tempDF = instAssessDataframe[['Longitude','Latitude','res']]
        saveModule.save_tiff_from_df(tempDF,['res'],0.0,path_out+str(instant)+'_anomalyInstantMap.tif',coordSystem)
    
#Building the dataFrame of "maps"
dictInstMaps['Longitude'] = instAssessDataframe.Longitude
dictInstMaps['Latitude'] = instAssessDataframe.Latitude
instMapsDataframe = pd.DataFrame(dictInstMaps)

#Regular expression for "cloud_\d" and "anomaly_\d"
listVars = instMapsDataframe.keys().to_list()
cloudVarList = [ 'cloud_' in i for i in listVars]
anomalyVarList = [ 'anomaly_' in i for i in listVars]

instMapsDataframe['freq'] = instMapsDataframe.loc[:,anomalyVarList].sum(axis=1) / instMapsDataframe.loc[:,cloudVarList].sum(axis=1)

#Save the output maps
if optionSaves[6]:
    tempDF = instMapsDataframe[['Longitude','Latitude','freq']]
    saveModule.save_tiff_from_df(tempDF,['freq'],0.0,path_out+str(instant)+'_anomalyMap.tif',coordSystem)
    
print('End of process...')




# instAssessDataframe = MFM_DataBase.build_assessment_instant_dataframe__withClouds_andRefModis(assessmentImage,
#                                                           dictImageInstants,
#                                                           cloudMasksInstants,
#                                                           int(instant),
#                                                           listAtts,
#                                                           geometry,
#                                                           scale,
#                                                           coordSystem,sensorName,baseModisRef)