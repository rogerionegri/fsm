#---------------------------------------------
#MÓDULO COM CONFIGURAÇÕES --------------------
#---------------------------------------------

def get_sensor_name():

    sensorName = {'multispec': 'landsat',
                  'multispecHR': 'sentinel',
                  'multispecLR': 'modis',
                  'precipitation': 'chirps',
                  'temperature': 'era5',
                  'fire': 'firms',
                  'humidity': 'fldas',          
                   } 

    return sensorName