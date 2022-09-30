# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:20:35 2022

@author: sergio.morales
"""

def getLastEv(debug=False):
    pathlastloc='runtime/lastloc.txt'
    import pandas as pd
    import ovdas_getfromdb_lib as gdb
    import datetime as dt
    fini = dt.datetime.strftime(dt.datetime.utcnow() - dt.timedelta(days=1), '%Y-%m-%d')
    ffin = dt.datetime.strftime(dt.datetime.utcnow() + dt.timedelta(days=1), '%Y-%m-%d')
    ev = pd.DataFrame(gdb.extraer_eventos(inicio=fini,final=ffin,volcan='*',ml='>0'))
    ev = ev[ev.tipoevento.isin(['VT','LP','LV','TO','HB'])]
    ev = ev[ev.cod!='-']
    ev = ev.tail(1)
    last_loc_hypo = str(ev['idevento'].iloc[0])+'-'+str(ev.fecha.iloc[0].year)
    print("Último evento localizado por Hypo: "+last_loc_hypo)
    import os
    if debug==False:
        if os.path.isfile(pathlastloc):
            last_loc_NLL = open(pathlastloc).read().splitlines()[0]   
            with open(pathlastloc, "w") as text_file:
                print(f"{last_loc_hypo}", file=text_file)
        else:
            print('Archivo lastloc.txt no existe, creando...')
            with open(pathlastloc, "w") as text_file:
                print("None", file=text_file)
            last_loc_NLL = open(pathlastloc).read().splitlines()[0] 
        if last_loc_hypo==last_loc_NLL:
            print('No hay eventos nuevos localizados')
            reloc=False
        else:
            print('A localizar!')    
            reloc=True
            with open(pathlastloc, "w") as text_file:
                print(f"{last_loc_hypo}", file=text_file)
    else:
        reloc=True
    return reloc, ev


def saveNLLcsv(voldf,ev,NLL,exec_point):
    if exec_point == 'local':
        raiz='//172.16.40.102/Monitoreo/ovv2023/csv/'
        
    elif exec_point == 'server':
        raiz='/mnt/puntocientodos/ovv2023/' 
    import pandas as pd
    import os
    volcan=voldf.nombre_db.iloc[0]
    fecha=str(ev.fecha.dt.year.iloc[0])+str(ev.fecha.dt.month.iloc[0]).zfill(2)+str(ev.fecha.dt.day.iloc[0]).zfill(2)
    #Guardar
    file = fecha+'_'+volcan+"_NLLoc.csv"
    row = (str(ev.fecha.iloc[0])+' '+str(NLL[0])+' '+str(NLL[1])+' '+str(NLL[2])+' '+str(NLL[3])+' '+str(NLL[4])+' '+str(NLL[5])
                                +' '+str(ev.tipoevento.iloc[0])+' '+str(ev.idevento.iloc[0])+'-'+str(ev.fecha.iloc[0].year)
                                )
    df = pd.DataFrame([row])
    raiz= raiz+str(ev.fecha.dt.year.iloc[0])+'/'+str(ev.fecha.dt.month.iloc[0]).zfill(2)+'/'+str(ev.fecha.dt.day.iloc[0]).zfill(2)+'/'
    local = 'csv/'+str(ev.fecha.dt.year.iloc[0])+'/'+str(ev.fecha.dt.month.iloc[0]).zfill(2)+'/'+str(ev.fecha.dt.day.iloc[0]).zfill(2)+'/'
    print(raiz)
    if (os.path.exists(raiz+file)==True):     
        df.to_csv(raiz+file,sep=';',index=False,mode='a',header=False)
    else: 
        if (os.path.exists(raiz)==True):
            df.to_csv(raiz+file,sep=';',index=False,header=False)
        else:
            os.makedirs(raiz)
            df.to_csv(raiz+file,sep=';',index=False,header=False)

    if (os.path.exists(local+file)==True):     
        df.to_csv(local+file,sep=';',index=False,mode='a',header=False)
    else: 
        if (os.path.exists(local)==True):
            os.makedirs(local)
            df.to_csv(local+file,sep=';',index=False,header=False)
        else:
            os.makedirs(local)
            df.to_csv(local+file,sep=';',index=False,header=False)