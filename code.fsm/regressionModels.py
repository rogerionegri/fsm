#Módulo para definição de modelos de regressão e afins
#=====================================================

#--------------------------------------------
#IMPORTS ------------------------------------
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RandomizedSearchCV
from sklearn.ensemble import IsolationForest
from sklearn import svm
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import make_scorer
import scipy


#--------------------------------------------
def get_model(modelName,nJobs):
    if modelName == 'RFR':
        random_grid = {'bootstrap': [True],                 
                       'max_depth': [10, 50, 100, None],    
                       'max_features': ['auto'],            
                       'n_estimators': [130, 180, 230]}
        model = RandomForestRegressor()
        model_RandSearh = RandomizedSearchCV(estimator = model, param_distributions = random_grid, n_iter = 20, cv = 3, verbose=2, random_state=0, n_jobs = nJobs, scoring=make_scorer(cohen_kappa_score))

    if modelName == 'SVM':
        # Grade de parâmetros - RandomizedSearchCV 
        random_grid = {'kernel': ['rbf'],
                       'gamma': scipy.stats.expon(scale=.100),
                       'C': [0.01,0.1,1,10,100,1000,10000]
                       }
        model = svm.SVR()
        model_RandSearh = RandomizedSearchCV(estimator = model, param_distributions = random_grid, n_iter = 20, cv = 3, verbose=2, random_state=0, n_jobs = nJobs, scoring=make_scorer(cohen_kappa_score))

    return model_RandSearh

#--------------------------------------------
def get_detector(modelName,nJobs):

    if modelName == 'OCSVM':
    
        random_grid = {'kernel': ['rbf'],
                       'gamma': scipy.stats.expon(scale=.100),
                       'nu': scipy.stats.expon(scale=.100),
                      }
        mmodel_adjustRSodel = svm.OneClassSVM(shrinking=True, kernel='rbf', nu=0.001)
       

    if modelName == 'IsolationForest':
        random_grid = {'bootstrap': [True],
                       'contamination': [0, 0.01, 0.05, 0.1],
                       'n_estimators': [100, 150, 200]
                      }
        model_adjustRS = IsolationForest(n_estimators=150,contamination=0.01)


    if modelName == 'svm':
    
        random_grid = {'kernel': ['rbf'],
                       'gamma': scipy.stats.expon(scale=.100),
                       'penalty': scipy.stats.expon(scale=-4),
                      }
        model_adjustRS = svm.SVC(shrinking=True, kernel='rbf', gamma=0.001, C=1000)


    return model_adjustRS



       