# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:37:29 2022

@author: sergio.morales
"""

def doNLL(evsel,exec_point,fini='',ffin='',signature='sergi'):
    idev = evsel.idevento.iloc[0]
    year = evsel.fecha.dt.year.iloc[0]    
    VpVs=1.78
    yNum,zNum=2000,1050
    ro = 2.7 #Densidad
    gridsize = 0.1 #Initial cell size in km
    import pandas as pd
    pd.options.mode.chained_assignment = None  # default='warn'
    import ovdas_getfromdb_lib as gdb
    import datetime as dt
    VpVs=1.78
    yNum,zNum=2000,1050
    ro = 2.7 #Densidad
    gridsize = 0.1 #Initial cell size in km
    
    from obspy import read_events
    hypo, NLL,voldf = mainNLL(exec_point,yNum,zNum,gridsize,VpVs,ro,evsel,signature)
    cat = read_events("loc/last.hyp") 
    lonN=cat[0].origins[0].longitude
    latN=cat[0].origins[0].latitude
    profN=cat[0].origins[0].depth/1000
    erlonN = cat[0].origins[0].longitude_errors.uncertainty
    erlatN = cat[0].origins[0].latitude_errors.uncertainty
    erzN=cat[0].origins[0].depth_errors.uncertainty/1000
    volnet = gdb.get_metadata_wws(volcan='*')
    volnet = volnet[volnet.tipo=='SISMOLOGICA']
    estas,fases =[],[]
    for i in range(0,len(cat[0].picks)):
        esta = cat[0].picks[i].waveform_id.station_code
        fase = cat[0].picks[i].phase_hint
        estas.append(esta[:-1])
        fases.append(fase)
    estadf = pd.DataFrame([estas,fases]).T
    NLL=[latN,lonN,profN,erlonN,erlatN,erzN]
    return NLL,estadf,hypo,voldf


def mainNLL(exec_point,yNum,zNum,gridsize,VpVs,ro,evsel,signature,verbose=True):
    import ovdas_getfromdb_lib as gdb
    import subprocess
    from random import randint
    volcanes =gdb.get_metadata_volcan('*',rep='y')
    volcanes = volcanes.drop_duplicates(subset='nombre', keep="first")
    print(evsel)
    voldf = volcanes[volcanes.id==evsel.idvolc.iloc[0]]
    voldf.latitud=voldf.latitud
    voldf.longitud=voldf.longitud
    if voldf.id.iloc[0]==40:
        voldf['altitud']=2500
        voldf['nref']=2500
    randomint = str(randint(1,100000))
    path1,lasttop = create_Vel2Grid_input(exec_point,voldf,yNum,zNum,gridsize,VpVs,ro,randomint) #CREAR Input para Vel2Grid
    import platform
    if platform.system()=='Windows':
        proc = subprocess.run(['wsl','Vel2Grid',path1], stdout=subprocess.PIPE, universal_newlines=True) #EJECUTA Vel2Grid  
        print(proc.stdout)
    else:
        proc = subprocess.run(['Vel2Grid',path1]) #EJECUTA Vel2Grid  
    path2,stacods = create_Grid2Time_input(evsel,voldf,randomint)#CREAR Input para Grid2Time
    if platform.system()=='Windows':
        proc = subprocess.run(['wsl','Grid2Time',path2], stdout=subprocess.PIPE, universal_newlines=True) # EJECUTA Grid2Time
        print(proc.stdout)
    else:
        subprocess.run(['Grid2Time',path2]) # EJECUTA Grid2Time
    create_obs_input(evsel) #Create observations file
    path3 = create_NLLoc_input(randomint,voldf,signature,evsel,stacods,VpVs) #Create Input para NLLloc
    if platform.system()=='Windows':
        proc = subprocess.run(['wsl','NLLoc',path3], stdout=subprocess.PIPE, universal_newlines=True) # EJECUTA NLL
        print(proc.stdout)
    else:
        subprocess.run(['NLLoc',path3]) # EJECUTA NLL        
    voldf.latitud=voldf.latitud
    voldf.longitud=voldf.longitud
    idev=evsel.idevento.iloc[0]
    if platform.system()=='Windows':    
        subprocess.run(['wsl','LocSum','./loc/'+str(idev)+'.sum.grid0.loc','1','./loc/'+str(idev),'./loc/'+str(idev)+'.*.*.grid0.loc']) # EJECUTA LocSum
    else:
        subprocess.run(['LocSum','./loc/'+str(idev)+'.sum.grid0.loc','1','./loc/'+str(idev),'./loc/'+str(idev)+'.*.*.grid0.loc']) # EJECUTA LocSum        
    if platform.system()=='Windows':   
        subprocess.run(['wsl','Grid2GMT',path1,'./loc/'+str(idev),'./gmt','L','S']) # EJECUTA LocSum
    else:
        subprocess.run(['Grid2GMT',path1,'./loc/'+str(idev),'./gmt','L','S']) # EJECUTA LocSum
    import os
    try:
        os.remove('gmt'+str(idev)+'.LS.gmt')
        os.remove('gmt'+str(idev)+'.LS.ps')
    except:
        ()
    evHypo71,latNLL,lonNLL,profNLL = get_loc_final(evsel,voldf) #Obtiene hipocentro
    
    
    
    #fig_NLL_single(evsel)
    if verbose==True:
        print('Hypo71    : Latitud: '+str(float(evHypo71.latitud.iloc[0]))+' , Longitud: '+str(float(evHypo71.longitud.iloc[0]))+', Profundidad: '+str(float(evHypo71.profundidad_abs.iloc[0]))+' km')
        print('NonLinLoc : Latitud: '+str(latNLL)+' , Longitud: '+str(lonNLL)+', Profundidad: '+str(profNLL)+' km')
    return evHypo71,[latNLL,lonNLL,float(profNLL)],voldf

blacklist=['NADA']
blacklist2=['NADA']

def delayf(altitud, NREF, volcan):
    import numpy as np
    #volcanes del norte excepto láscar
    if volcan in (1,2,3,4,5,6,7,8,9,11):
        mod_1D_v = [4.00,4.80,5.95,6.45,6.61,6.91,7.10,7.53,7.71,8.11]
        mod_1D_z = [0,2.3,4.7,14.7,24.7,34.7,44.7,54.7,64.7,104.7]
    #Láscar
    elif volcan == 10:
        mod_1D_v = [1.00,1.30,1.80,2.20,2.70,3.00,3.50,4.00,4.50,5.00,6.00,6.50,7.00]
        mod_1D_z = [0,0.5,1.0,1.5,2,3,4,5,6,8,10,20,40]
    #Maule
    elif volcan==19:
        mod_1D_v = [4.31,4.39,4.46,5.65,6.00,6.89,7.40,7.76,7.94,8.34]
        mod_1D_z = [0,4.0,5.0,8.0,18,23,38,48,58,93]
    #Chillán
    elif volcan==21:
        mod_1D_v = [3.60,3.80,4.00,4.20,4.40,4.80,5.20,5.60,6.28,6.89,7.40,7.76,7.94,8.34]
        mod_1D_z = [0,0.5,1.0,1.5,2,3,4,6,8,23,38.8,48,58,93]
    #Copahue
    elif volcan==23:
        mod_1D_v = [4.00,4.50,4.80,5.30,5.70,6.00,6.30,6.89,7.40,7.76,7.94,8.34]
        mod_1D_z = [0,0.9,2.9,5.4,7.9,12.9,17.9,22.9,37.9,47.9,57.9,92.9]
    #Osorno
    elif volcan==36:
        mod_1D_v = [3.90,4.00,4.20,4.40,4.60,4.80,5.00,5.20,5.40,5.60,5.80,6.00,6.30,6.50,6.70,6.89,7.40,7.76,7.94,8.34]
        mod_1D_z = [0,0.5,1.0,1.5,2,2.5,3,3.5,4,5,6,7,8,10,15,20,35,45,55,90]
    #Todo el resto
    else:
        mod_1D_v = [4.39,5.51,6.28,6.89,7.40,7.76,7.94,8.34]
        mod_1D_z = [0,2,5,20,35,45,55,90]
    ######################################PROCESAMIENTO
    h_esta = (NREF-altitud)/1000
    i=0
    for n in range(0,len(mod_1D_z)):
        if (h_esta>mod_1D_z[n]) and (h_esta<=mod_1D_z[n+1]):
            k = i
        i=i+1
    if h_esta>0:
        delay = 0
        if k>0: #EXISTEN CAPAS COMPLETAS
            #SUMADO DE CAPAS COMPLETAS
            for x in range(0,k):
                grosor =  (mod_1D_z[x+1]-mod_1D_z[x])
                vel = mod_1D_v[x]
                delay_fullcapa = grosor/vel
                delay=delay+delay_fullcapa
        else:
            x=-1
        resto_h = h_esta-mod_1D_z[x+1]
        vel_resto = mod_1D_v[x+1]
        delay_resto=resto_h/vel_resto
        #SUMADO DE CAPA INCOMPLETA
        delay_total = delay+delay_resto
        d=np.round(delay_total,2)*-1
    else:
        grosor = -h_esta
        vel = mod_1D_v[0]
        delay = grosor/vel
        d = np.round(delay,2)*1
    return d

def get_delays(stacods,voldf,vPvS):
    import sys
    sys.path.append('//172.16.40.10/sismologia/pyovdas_lib/')
    import ovdas_getfromdb_lib as gdb
    estadf = gdb.get_metadata_wws(volcan='*')
    delays=[]
    for est in stacods:
        est = est[1:]
        try:
            altitud = estadf[estadf.codcorto==est].altitud.iloc[0]
            NREF = voldf.nref.iloc[0]
            delayP=delayf(altitud, NREF,voldf.id.iloc[0])
            delayS=delayP/vPvS
            delays.append('LOCDELAY '+est+'Z P 1 '+str(delayP))
            delays.append('LOCDELAY '+est+'Z S 1 '+str(delayS))
        except:
            print('Estación '+est+' no existe en db')
    return delays


def create_common_variables(voldf,randomint):
    com_CONTROL = "CONTROL 1 "+randomint   
    com_TRANS = "TRANS AZIMUTHAL_EQUIDIST WGS-84 "+str(float(voldf.latitud.iloc[0]))+' '+str(float(voldf.longitud.iloc[0]))+' 0'
    return com_CONTROL,com_TRANS

def get_layers(exec_point,voldf,VpVs,ro):
    import numpy as np
    if exec_point=='server':
        root = '/mnt/puntodiez/pyovdas_lib/mod/'
    elif exec_point=='local':
        root = '//172.16.40.10/sismologia/pyovdas_lib/mod/'
    lines = open(root+str(voldf.id.iloc[0])+'.txt').read().splitlines()[:-2]
    layers = []
    for line in lines:
        line = line.split(' ')
        prof = float(line[-1])-voldf.nref.iloc[0]/1000
        linea='LAYER '+str(prof)+' '+line[2]+ ' 0.00 '+str(np.round(float(line[2])/VpVs,2))+' 0.0 '+str(ro)+ ' 0.0'
        layers.append(linea)
    return layers,line[-1]

def create_Vel2Grid_input(exec_point,voldf,yNum,zNum,gridsize,VpVs,ro,randomint):  
    outdir_VEL2Grid_input= './inputs/1_Vel2Grid_input_'+voldf.nombre_db.iloc[0]
    com_CONTROL, com_TRANS = create_common_variables(voldf,randomint)
    com_VGOUT   = "VGOUT ./model/layer"
    com_VGTYPE  = "VGTYPE P"
    volpeak = float(voldf.altitud.iloc[0]/1000)*-1
    com_VGGRID  = ("VGGRID 2 "+str(yNum)+' '+str(zNum)+' 0.0 0.0 '+str(volpeak)+' ' 
                   +str(gridsize)+' '+str(gridsize)+' '+str(gridsize)+' SLOW_LEN')
    com_LAYERS,lasttop = get_layers(exec_point,voldf,VpVs,ro)
    file = open(outdir_VEL2Grid_input, "w") 
    with open(outdir_VEL2Grid_input, 'a') as the_file:
        the_file.write(com_CONTROL+'\n'+com_TRANS+'\n'+com_VGOUT+'\n'+com_VGTYPE+'\n'+com_VGGRID)
        for line in com_LAYERS:
            the_file.write('\n'+line)
        file.close()
    return outdir_VEL2Grid_input,lasttop

def get_stas(evsel):
    import ovdas_getfromdb_lib as gdb
    estadf = gdb.get_metadata_wws(volcan='*',estadoInst='*',estadoEst='*')
    fases = gdb.get_fasesLoc(evsel.idevento.iloc[0], evsel.fecha.iloc[0].year)
    fases = fases[~fases.cod.isin(blacklist2)]

    print(fases)
    staline = []
    stacods =[]

    for index,row in fases.iterrows():
        stacod = str(row.cod)
        try:
            stalat = str(estadf[estadf.codcorto==stacod[1:]].latitud.iloc[0])
            stalon = str(estadf[estadf.codcorto==stacod[1:]].longitud.iloc[0])
            stahei = str(estadf[estadf.codcorto==stacod[1:]].altitud.iloc[0]/1000)
            stacods.append(stacod)
        except:
            print('problemas, estación '+stacod+' no está en metadata_wws')
        staline.append('GTSRCE '+stacod[1:]+'Z LATLON '+stalat+' '+stalon+' 0.0 '+stahei)
    return staline,stacods

def create_Grid2Time_input(evsel,voldf,randomint):
    outdir_Grid2Time_input= './inputs/2_Grid2Time_input_'+voldf.nombre_db.iloc[0]
    com_CONTROL, com_TRANS = create_common_variables(voldf,randomint)   
    com_GTFILES = 'GTFILES  ./model/layer  ./time/layer P'
    com_GTMODE =  'GTMODE GRID2D ANGLES_YES'
    com_GTSRCE,stacods = get_stas(evsel)
    com_GT_PLFD = 'GT_PLFD  1.0e-3  0'
    file = open(outdir_Grid2Time_input, "w") 
    with open(outdir_Grid2Time_input, 'a') as the_file:
        the_file.write(com_CONTROL+'\n'+com_TRANS+'\n'+com_GTFILES+'\n'+com_GTMODE)
        for line in com_GTSRCE:
            the_file.write('\n'+line)
        the_file.write('\n'+com_GT_PLFD)
        file.close()   
    return outdir_Grid2Time_input,stacods
    
def create_obs_input(evsel):
    outdir_obs_input = './obs/'+str(evsel.idevento.iloc[0])+'.obs'
    import sys
    from datetime import datetime
    #import pandas as pd
    sys.path.append('//172.16.40.10/sismologia/pyovdas_lib/')
    import ovdas_getfromdb_lib as gdb
    fases = gdb.get_fasesLoc(evsel.idevento.iloc[0], evsel.fecha.iloc[0].year)
    fases = fases[~fases.cod.isin(blacklist2)]
    file = open(outdir_obs_input, "w") 
    with open(outdir_obs_input, 'a') as the_file:

        for index,row in fases.iterrows():
            if row.fullstring is None:
               if row.pesop==0:pesop=0.01
               elif row.pesop==1:pesop=0.5
               elif row.pesop==2:pesop=1.0
               elif row.pesop==3:pesop=2.0
               elif row.pesop==4:pesop=3.0
               elif row.pesop==9:pesop=999999.0
               
               if row.pesos==0:pesos=0.01
               elif row.pesos==1:pesos=0.5
               elif row.pesos==2:pesos=1.0
               elif row.pesos==3:pesos=2.0
               elif row.pesos==4:pesos=3.0
               elif row.pesos==9:pesos=999999.0 
               
               polaridad=row.polaridad
               if polaridad==" ":polaridad='?'
               tiempop=datetime.fromtimestamp(row.tiempop)
               fila = (row.cod[1:]+'Z ? ? '+row.tipoarribo.lower()+' P '+polaridad+' '+str(tiempop.year)+str(tiempop.month).zfill(2)+
                       str(tiempop.day).zfill(2)+' '+str(tiempop.hour).zfill(2)+str(tiempop.minute).zfill(2)+' '+
                       str(tiempop.second).zfill(2)+'.'+str(int(tiempop.microsecond/10000))+ 
                       " GAU "+
                       str(pesop)+
                       " -1.00e+00 -1.00e+00 -1.00e+00"
                       )
               the_file.write(fila+'\n')
               #if evsel.idvolc.iloc[0] != 40:
               if row.tiempos>0:
                    tiempos=datetime.fromtimestamp(row.tiempos)
                    fila = (row.cod[1:]+'Z ? ? '+row.tipoarribo.lower()+' S ? '+str(tiempos.year)+str(tiempos.month).zfill(2)+
                            str(tiempos.day).zfill(2)+' '+str(tiempos.hour).zfill(2)+str(tiempos.minute).zfill(2)+' '+
                            str(tiempos.second).zfill(2)+'.'+str(int(tiempos.microsecond/10000))+ 
                            " GAU "+
                       str(pesos)+
                            " -1.00e+00 -1.00e+00 -1.00e+00"
                            )                   
                    the_file.write(fila+'\n')
            else:
                polaridad=row.fullstring[6]
                if polaridad==" ":polaridad='?'
                fila = (row.fullstring[0:4]+" ? ? "+row.fullstring[4].lower()+" "+row.fullstring[5]+" "+polaridad+
                      ' 20'+row.fullstring[9:15]+" "+row.fullstring[15:19]+" "+row.fullstring[19:24]+
                      " GAU "+
                       str(pesop)+
                      " -1.00e+00 -1.00e+00 -1.00e+00")
                the_file.write(fila+'\n')
                if (row.fullstring[31:36]) !='00.00':
                    filas = (row.fullstring[0:4]+" ? ? "+row.fullstring[4].lower()+" "+row.fullstring[37]+" ? 20"+
                          row.fullstring[9:15]+" "+row.fullstring[15:19]+" "+row.fullstring[31:36]+
                          " GAU" +
                       str(pesos)+
                          " -1.00e+00 -1.00e+00 -1.00e+00")
                    the_file.write(filas+'\n')
            
                    the_file.write(row.fullstring+'\n')
        file.close()   

def create_NLLoc_input(randomint,voldf,signature,evsel,stacod,vPvS):
    outdir_NLLoc_input = './inputs/4_NLLoc_input_'+voldf.nombre_db.iloc[0]
    com_CONTROL, com_TRANS = create_common_variables(voldf,randomint)  
    com_LOCSIG      = 'LOCSIG '+str(signature)
    com_LOCCOM      = 'LOCCOM Localizacion de evento id:'+str(evsel.idevento.iloc[0])
    com_LOCFILES    = 'LOCFILES ./obs/'+str(evsel.idevento.iloc[0])+'.obs NLLOC_OBS ./time/layer  ./loc/'+str(evsel.idevento.iloc[0])
    com_LOCHYPOUT   = 'LOCHYPOUT SAVE_NLLOC_ALL  SAVE_HYPOINV_SUM'
    com_LOCSEARCH   = 'LOCSEARCH  OCT 10 10 10 0.001 50000 5000 0 1'
    volpeak = float(voldf.altitud.iloc[0]/1000)*-1
    com_LOCGRID     = 'LOCGRID  101 101 101  -50.0 -50.0 '+str(volpeak)+'  1 1 1   PROB_DENSITY  SAVE'
    com_LOCMETH     = 'LOCMETH EDT_OT_WT 9999.0 4 -1 -1 1.68 6 -1.0 1'
    com_LOCGAU      = 'LOCGAU 0.2 0.0'
    com_LOCGAU2     = 'LOCGAU2 0.01 0.05 0.9'
    com_PHASEID1    = 'LOCPHASEID  P   P p G PN PG'
    com_PHASEID2    = 'LOCPHASEID  S   S s G SN SG'
    com_LOCANGLES   = 'LOCANGLES ANGLES_YES 5'
    com_LOCQUAL2ERR = 'LOCQUAL2ERR 0.1 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 99999.9'
    com_LOCMAG      = 'LOCMAG ML_HB 1.0 1.110 0.00189'
    com_LOCPHSTAT   = 'LOCPHSTAT 9999.0 -1 9999.0 1.0 1.0 9999.9 -9999.9 9999.9'
    com_LOCSTAWT    = 'LOCSTAWT 1 -1'
    #com_LOCTOPO_SURFACE = 'LOCTOPO_SURFACE '+'./grd/'+str(voldf.cod.values[0])+'.nc.asc'+' 0'
    #com_LOCDELAY = get_delays(stacod,voldf,vPvS)
    file = open(outdir_NLLoc_input, "w") 
    with open(outdir_NLLoc_input, 'a') as the_file:
        the_file.write(com_CONTROL+'\n'+com_TRANS+'\n'+com_LOCSIG+'\n'+com_LOCCOM+'\n'+com_LOCFILES+'\n'+com_LOCHYPOUT+'\n'+com_LOCSEARCH+'\n'+com_LOCGRID+'\n'+
                       com_LOCMETH+'\n'+com_LOCGAU+'\n'+com_LOCGAU2+'\n'+com_PHASEID1+'\n'+com_PHASEID2+'\n'+com_LOCANGLES+'\n'+
                       #com_LOCMAG+'\n'+com_LOCPHSTAT+'\n'+com_LOCTOPO_SURFACE)
                       com_LOCMAG+'\n'+com_LOCPHSTAT+'\n'+com_LOCQUAL2ERR+'\n'+com_LOCSTAWT)
        #for line in com_LOCDELAY:
        #    the_file.write('\n'+line)
        file.close()  
    return outdir_NLLoc_input

def get_loc_final(evsel,voldf):
    evento=evsel
    file =(str(evento.idevento.iloc[0])+'.sum.grid0.loc.hyp')
    #file =(str(evento.idevento.iloc[0])+'.19'+str(evento.fecha.iloc[0])[2:4]+str(evento.fecha.iloc[0])[5:7]+
    #           str(evento.fecha.iloc[0])[8:10]+'.'+str(evento.fecha.iloc[0])[11:13]+str(evento.fecha.iloc[0])[14:16]+
    #           str(evento.fecha.iloc[0])[17:19]+'.grid0.loc.hyp')
    #file='last.hyp'
    with open("./loc/"+file) as file_in:
        lines = []
        for line in file_in:
            lines.append(line)
        latNLL = float((lines[6].split(' ')[12])) 
        lonNLL = float((lines[6].split(' ')[14])) 
        profNLL= (lines[6].split(' ')[16])[:-2] 
    evsel['profundidad_abs']=evsel.profundidad.iloc[0]-voldf.nref.iloc[0]/1000
    return evsel,latNLL,lonNLL,float(profNLL)