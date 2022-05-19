#--------------------------------------------
#IMPORTS ------------------------------------
import numpy as np
from sklearn.mixture import GaussianMixture
import plotAux


def gmm_n_components_selection(dataset):
    #Adaptado de https://stackoverflow.com/questions/39920862/model-selection-for-gaussianmixture-by-using-gridsearch

    bic = np.zeros(10)
    n = np.arange(1,11)
    models = []
    #loop through each number of Gaussians and compute the BIC, and save the model
    for i,j in enumerate(n):
        #create mixture model with j components
        gmm = GaussianMixture(n_components=j, random_state=0)
        #fit it to the data
        gmm.fit(dataset)
        #compute the BIC for this model
        bic[i] = gmm.bic(dataset)
        #add the best-fit model with j components to the list of models
        models.append(gmm)

    #pos = np.argmin(bic)
    pos = np.argmax( (np.roll(bic,1) - bic)/np.abs(bic) )
    #pos = np.where( (np.roll(bic,1) - bic)/np.abs(bic) > 10**(-3) )

    return bic[pos], models[pos]

#--------------------------------------------
def fire_spots_df(tabRaw,ind):#,gmmComponents):

    #Criando uma cópia, já que a atribuição de dataFrames é feita por referência!
    tabFire = tabRaw.copy()

    #Selecao baseada no "nome da coluna"
    coi = [col for col in tabFire.columns if 'Foco' in col]  #...column of interest
    tabFire['SomaFocos'] = tabFire.loc[:,coi].sum(axis=1)
    
    #Seleciona as linhas com valor não nulo em "SomaFocos"
    tabFireOccurence = tabFire.loc[tabFire['SomaFocos'] > 0]

    
    if tabFireOccurence.shape[0] > 0:
    
        #Descarta as colunas denecessárias deste dataFrame (que é apenas intermediário...)
        tabFireOccurence.drop( tabFireOccurence.loc[:,coi].columns , axis=1, inplace=True )#.reset_index(drop=True)

        #Buscar repetições e replicar a mesma quantidade de vezes
        tabFireOccurenceRep = tabFireOccurence.loc[tabFireOccurence.index.repeat(tabFireOccurence.SomaFocos)].reset_index(drop=True)

        #Cada linha da tabela original "tabFire" será comparada com os elementos de "tabFireOccurenceRep"
        tabFire['logProb'] = 0.0
        tabFire['probAprox'] = 0.0
        tabFire['prob'] = 0.0

        #definir o "dataset"
        dataset = (tabFireOccurenceRep[['Latitude','Longitude']]).to_numpy()
        
        #Determinar a distribuição de probabilidade da ocorência de fogo via GMM
        #Em "gmm_n_components_selection" é feita a escolha do número de componentes e ajuste do modelo
        #...está fixo aqui que a busca vai até 10 componentes  
        bic, gm = gmm_n_components_selection(dataset)  
        
        #Conjunto de "pontos" completo
        predictSet = (tabFire[['Latitude','Longitude']]).to_numpy()
    
        #proba = np.exp( kde.score_samples( predictSet ) )
        proba = np.exp( gm.score_samples( predictSet ) )   #A saída de .score_samples() é uma "Log-probabilidade"
        
        deltaX = predictSet[:,0].max() - predictSet[:,0].min()
        deltaY = predictSet[:,1].max() - predictSet[:,1].min()
        cellArea = ( deltaX/np.unique( predictSet[:,0] ).shape[0] ) * ( deltaY/np.unique( predictSet[:,1] ).shape[0] )

        #Vou deixar as duas opções
        tabFire['logProb'] = proba
        tabFire['probAprox'] = proba*cellArea #É uma aproximação grosseira, pois (cellArea * proba).sum() é diferente de 1
        tabFire['prob'] = proba/proba.max() #É uma aproximação grosseira, pois (cellArea * proba).sum() é diferente de 1

    else:
        #Caso onde não são observados focos de incêndio
        tabFire['logProb'] = 0.0
        tabFire['probAprox'] = 0.0
        tabFire['prob'] = 0.0


    #Temporario...
    plotAux.plot_check_map(tabFire,'/home/rogerio/GIT/msc.andrea/Codes/roger/saidasTeste/prob'+ind+'.png','prob')
    print('ploted prob')
    plotAux.plot_check_map(tabFire,'/home/rogerio/GIT/msc.andrea/Codes/roger/saidasTeste/logProb'+ind+'.png','logProb')
    print('ploted log prob')

    return tabFire

#--------------------------------------------
def fire_spots_df_v2(tabRaw,ind):

    #Criando uma cópia, já que a atribuição de dataFrames é feita por referência!
    tabFire = tabRaw.copy()
    
    tabFire['indicador'] = 0 #Cria a coluna do rotulo de "ocorrência"

    #Selecao baseada no "nome da coluna"
    coi = [col for col in tabFire.columns if 'Foco' in col]  #...column of interest
    tabFire['SomaFocos'] = tabFire.loc[:,coi].sum(axis=1)
    
    #Seleciona as linhas com valor não nulo em "SomaFocos"
    tabFireOccurence = tabFire.loc[tabFire['SomaFocos'] > 0]    

    #   retornar uma tabela vazia mesmo!
    if tabFireOccurence.shape[0] > 0:

        #Buscar repetições e replicar a mesma quantidade de vezes
        tabFireOccurence = tabFireOccurence.loc[tabFireOccurence.index.repeat(tabFireOccurence.SomaFocos)].reset_index(drop=True)
        tabFireOccurence['indicador'] = 1
    
    #Temporario...
    plotAux.plot_check_map(tabFireOccurence,'/home/rogerio/GIT/msc.andrea/Codes/roger/saidasTeste/indica_'+ind+'.png','indicador')
    print('ploted indicador...')
    
    return tabFireOccurence

#--------------------------------------------
#Versão baseada nos valores (não normalziados) de NBR e NDWI
def fire_spots_df_v3(tabRaw,ind):

    #Criando uma cópia, já que a atribuição de dataFrames é feita por referência!
    tabFire = tabRaw.copy()
    
    tabFire['indicador'] = 0 #Cria a coluna do rotulo de "ocorrência"

    #Selecao baseada no "nome da coluna"
    coi = [col for col in tabFire.columns if 'Foco' in col]  #...column of interest
    tabFire['SomaFocos'] = tabFire.loc[:,coi].sum(axis=1)
    
    #Seleciona as linhas com valor não nulo em "SomaFocos"
    tabFireOccurence = tabFire.loc[tabFire['SomaFocos'] > 0]    

    #   retornar uma tabela vazia mesmo!
    if tabFireOccurence.shape[0] > 0:

        #Buscar repetições e replicar a mesma quantidade de vezes
        tabFireOccurence = tabFireOccurence.loc[tabFireOccurence.index.repeat(tabFireOccurence.SomaFocos)].reset_index(drop=True)
        
        coiNBR = [col for col in tabFireOccurence.columns if 'NBR' in col] #...column of interest
        
        tabFireOccurence['deltaBurn'] = tabFireOccurence.loc[:,coiNBR].max(axis=1) - tabFireOccurence.loc[:,coiNBR].min(axis=1)
        
        #http://www.dsr.inpe.br/sbsr2015/files/p0104.pdf  << limiares do deltaNBR...
        condition = tabFireOccurence.deltaBurn > 0.2 #np.logical_and( tabFireOccurence.deltaNBR < -0.5 , tabFireOccurence.maxNDWI < 0.3 )
        
        tabFireOccurence['indicador'][condition] = 1
        
        
    plotAux.plot_check_map(tabFireOccurence,'/home/rogerio/GIT/msc.andrea/Codes/roger/saidasTeste/indica_'+ind+'.png','indicador')
    print('ploted indicador...')

    plotAux.boxplot_attributes(tabFireOccurence,'/home/rogerio/GIT/msc.andrea/Codes/roger/saidasTeste/boxplot_'+ind+'.png',['NDVI1','NDWI1','Temp1'])

    return tabFireOccurence      


#Versão baseada nos valores (não normalziados) de NBR e NDWI
def fire_spots_df_v4(tabRaw,ind):

    #Criando uma cópia, já que a atribuição de dataFrames é feita por referência!
    tabFire = tabRaw.copy()
    
    tabFire['indicador'] = 0 #Cria a coluna do rotulo de "ocorrência"

    #Selecao baseada no "nome da coluna"
    coi = [col for col in tabFire.columns if 'Foco' in col]  
    tabFire['SomaFocos'] = tabFire.loc[:,coi].sum(axis=1)
    
    #Seleciona as linhas com valor não nulo em "SomaFocos"
    tabFireOccurence = tabFire.loc[tabFire['SomaFocos'] > 0]    

    #Se a tabela de focos é vazia... paciência...
    #   retornar uma tabela vazia mesmo!
    if tabFireOccurence.shape[0] > 0:

        #Buscar repetições e replicar a mesma quantidade de vezes
        tabFireOccurence = tabFireOccurence.loc[tabFireOccurence.index.repeat(tabFireOccurence.SomaFocos)].reset_index(drop=True)

        condition = tabFireOccurence.deltaBurn > 0.1 #np.logical_and( tabFireOccurence.deltaNBR < -0.5 , tabFireOccurence.maxNDWI < 0.3 )
        
        tabFireOccurence['indicador'][condition] = 1   
    
    plotAux.plot_check_map(tabFireOccurence,'/home/rogerio/GIT/msc.andrea/Codes/roger/saidasTeste/indica_'+ind+'.png','indicador')
    print('ploted indicador...')
    
    plotAux.boxplot_attributes(tabFireOccurence,'/home/rogerio/GIT/msc.andrea/Codes/roger/saidasTeste/boxplot_'+ind+'.png',['NDVI1','NDWI1','Temp1'])

    return tabFireOccurence
    





