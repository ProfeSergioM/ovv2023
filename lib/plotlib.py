# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 16:08:02 2022

@author: sergio.morales
"""

blacklist=['NADA']
blacklist2=['NADA']


def obtener_evento(year,idev):
    import pandas as pd
    import ovdas_getfromdb_lib as gdb
    pd.options.mode.chained_assignment = None  # default='warn'
    eventos = gdb.extraer_eventos(str(year)+'-01-01',str(year)+'-12-31',volcan='*',ml='>0')
    eventos = pd.DataFrame(eventos)
    evsel = eventos[eventos.idevento==idev]
    volcanes =gdb.get_metadata_volcan('*',rep='y')
    volcanes = volcanes.drop_duplicates(subset='nombre', keep="first")
    voldf = volcanes[volcanes.id==evsel.idvolc.iloc[0]]  
    import xarray as xr
    import platform
    if platform.system()=='Windows':
        data = xr.load_dataset('C:/nc/'+str(voldf.id_zona.values[0])+'/'+str(voldf.cod.values[0])+'.nc')
    else:
        data = xr.load_dataset('/home/sismologia/nc/'+str(voldf.id_zona.values[0])+'/'+str(voldf.cod.values[0])+'.nc')
    dem = data.fillna(0)
    return evsel,voldf,dem

def plotNLLinfo(evsel,exec_point):
    idev = evsel.idevento.iloc[0]
    year = evsel.fecha.dt.year.iloc[0]    
    evsel,voldf,dem = obtener_evento(year, idev) #Obtener evento dede la DB
    evento=evsel
    file =(str(evento.idevento.iloc[0])+'.sum.grid0.loc.hyp')
    with open("./loc/"+file) as file_in:
        lines = []
        for line in file_in:
            lines.append(line)
        latNLL = float((lines[6].split(' ')[12])) 
        lonNLL = float((lines[6].split(' ')[14])) 
        profNLL= float((lines[6].split(' ')[16])[:-2] )
    evsel['profundidad_abs']=evsel.profundidad.iloc[0]-voldf.nref.iloc[0]/1000
    evHypo71 = evsel
    hypo=([evHypo71.longitud.iloc[0],evHypo71.latitud.iloc[0],evHypo71.profundidad_abs.iloc[0],
             evHypo71.erh.iloc[0],evHypo71.erz.iloc[0]])
    NLL=[lonNLL,latNLL,profNLL]
    plot_NLL(voldf,hypo,NLL,evsel,str(evento.idevento.iloc[0]),exec_point)
    
    
def plot_NLL(voldf,hypo,NLL,evsel,strid,exec_point):
    import numpy as np
    filename='loc/'+strid+'.hyp'
    with open(filename) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
    scatter=[]
    go=False
    for line in lines:
        if line[0:7]=='SCATTER':
            go=True
        if go==True:
            scatter.append(line)
    scatter2=scatter[1:-3]
    PDF=[]
    for line in scatter2:
        PDF.append(float(line.split()[3]))
    PDF = np.array(PDF) 
    fasedists=[]
    go=False
    for line in lines:
        if line[0:5]=='PHASE':
            go=True
        if line[0:9]=='END_PHASE':
            go=False
        if go==True:
            fasedists.append(line)
    sta,dis=[],[]
    
    for line in fasedists[1:]:
        sta.append(line.split()[0])
        dis.append(line.split()[22])
    import pandas as pd
    dDICT = dict(zip(sta, dis))

    from datetime import datetime
    for line in lines:
        if line[0:10]=='GEOGRAPHIC': 
            linea = line.split()
            ot = datetime(int(linea[2]),int(linea[3]),int(linea[4]),int(linea[5]),
                          int(linea[6]),int(linea[7][0:2]),int(linea[7][3:])
                          )
    idev = evsel.idevento.iloc[0]
    import xarray as xr
    import numpy as np
    import matplotlib
    from matplotlib.colors import LightSource
    import platform
    if platform.system()=='Windows':
        data = xr.load_dataset('C:/nc/'+str(voldf.id_zona.values[0])+'/'+str(voldf.cod.values[0])+'.nc')
    else:
        data = xr.load_dataset('/home/sismologia/nc/'+str(voldf.id_zona.values[0])+'/'+str(voldf.cod.values[0])+'.nc')
    data = data.fillna(0)
    import sys
    sys.path.append('/mnt/puntodiez/pyovdas_lib/')
    import ovdas_getfromdb_lib as gdb
    from obspy import read_events
    import geopy
    import geopy.distance
    cat = read_events("loc/last.hyp")
    erlonN =cat[0].origins[0].longitude_errors.uncertainty
    erlatN =cat[0].origins[0].latitude_errors.uncertainty
    erproN = cat[0].origins[0].depth_errors.uncertainty/1000
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        import matplotlib
        import numpy as np
        new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap
    import pandas as pd
    from math import cos, radians
    PDF_lonlatXY= pd.read_csv('loc/'+str(idev)+'.scat.lonlat.XY',header=None,sep=' ')
    PDF_lonlatXY = PDF_lonlatXY.rename(columns={0:'lon',1:'lat'})
    PDF_lonlatXZ= pd.read_csv('loc/'+str(idev)+'.scat.lonlat.XZ',compression=None,header=None,sep=' ')
    PDF_lonlatXZ = PDF_lonlatXZ.rename(columns={0:'lon',1:'prof'})
    df_PDF = pd.DataFrame([PDF_lonlatXZ['lon'],PDF_lonlatXY['lat'],PDF_lonlatXZ['prof']]).T
    df_PDF['PDF']=(PDF-min(PDF))/(max(PDF)-min(PDF))
    londelta = 0.2
    loncenter=(hypo[0]+NLL[0])/2
    latcenter=(hypo[1]+NLL[1])/2
    aspect_ratio = 1/cos(radians(latcenter))
    latdelta = londelta/aspect_ratio
    
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import cartopy.crs as ccrs
    fig = plt.figure(figsize=(16,9))
    
    grilla_externa = fig.add_gridspec(1,2, wspace=0.0, hspace=0.0)
    grilla_externa.update(left=0.05,right=0.9,top=0.85,bottom=0.05,wspace=0.3,hspace=0.09)
    mapaax = fig.add_subplot(grilla_externa[0])
    mapaax.axis('off')
    
    mapaloc = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec=mapaax, wspace=0.0, hspace=0.0,width_ratios=[5,9],height_ratios=[9,5])
    mapa1 = fig.add_subplot(mapaloc[1],projection=ccrs.PlateCarree())
    mapa1.set_title('Localización Hypo71/NonLinLoc (NLL)',fontsize=12)
    mapa1.set_extent([loncenter-londelta,loncenter+londelta,
                     latcenter-latdelta,latcenter+latdelta])    

    mapa1_proflat = fig.add_subplot(mapaloc[0])
    mapa1_proflat.set_ylim(mapa1.get_ylim())
    mapa1_proflat.grid(True,linestyle='--',alpha=0.5,linewidth=0.5,color='black')
    
    mapa1_proflon = fig.add_subplot(mapaloc[3])
    mapa1_proflon.set_xlim(mapa1.get_xlim())
    mapa1_proflon.grid(True,linestyle='--',alpha=0.5,linewidth=0.5,color='black')
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    nref=voldf.nref.iloc[0]
    zlim = 25
    import matplotlib.ticker as ticker
    hei = np.array(data['z'])
    cmap = matplotlib.cm.get_cmap('gist_earth')
    ls = LightSource(azdeg = 180, altdeg = 60)
    cmap = truncate_colormap(cmap)
    rgb = ls.shade(hei, cmap=cmap)
    mapa1.imshow(rgb,extent=(float(data['x'].min()),float(data['x'].max()),
                       float(data['y'].min()),float(data['y'].max())),
          interpolation='None', alpha=0.5,origin='lower')
  
    mapa1.set_aspect('auto')
    mapa1.set_aspect('auto')
    gl = mapa1.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=False)
    gl.xlocator = ticker.MultipleLocator(0.03)
    gl.ylocator = ticker.MultipleLocator(0.02)
    mapa1_proflon.set_ylim(zlim,-nref/1000)
    mapa1_proflat.set_xlim(zlim,-nref/1000)
    origin = geopy.Point(latitude=hypo[1],longitude=hypo[0])
    d_x = geopy.distance.distance(kilometers=hypo[3]*2)
    step = d_x.destination(point=origin, bearing=90)
    erlonH=step[1]-hypo[0]
    d_y = geopy.distance.distance(kilometers=hypo[4]*2)
    step = d_y.destination(point=origin, bearing=0)
    erlatH=step[0]-hypo[1]
    scat = mapa1.scatter(df_PDF['lon'],df_PDF['lat'],s=1,c=df_PDF['PDF'])    
    mapa1_proflat.scatter(df_PDF['prof'],df_PDF['lat'],s=1,c=df_PDF['PDF'])
    mapa1_proflon.scatter(df_PDF['lon'],df_PDF['prof'],s=1,c=df_PDF['PDF'])    
    
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes


    
    axins = inset_axes(mapa1,                  width="5%",  # width = 5% of parent_bbox width
                   height="50%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.0, 0.0, 1, 1),
                   bbox_transform=mapa1.transAxes,
                   borderpad=0
                   )
    
    cbar = fig.colorbar(scat, cax=axins, orientation="vertical")    
    cbar.set_label('Calidad de localización normalizada')
    #mapa1.plot(df_PDF.lon,df_PDF.lat,'.',ms=0.5,alpha=0.5,label='PDF')
    mapa1.plot(hypo[0],hypo[1],'*',mfc='C3',markersize=10,mec='k',zorder=999,label='Hypo71')
    mapa1.errorbar(hypo[0],hypo[1],yerr=erlatH,xerr=erlonH,
                   marker='*',color='C3',ecolor='C3',ms=0,alpha=1,lw=0,elinewidth=1,zorder=998)
    mapa1.plot(NLL[0],NLL[1],marker='*',mfc='C2',markersize=10,mec='k',lw=0,zorder=999,label='NonLinLoc')
    mapa1.errorbar(NLL[0],NLL[1],yerr=erlatN,xerr=erlonN,
                   marker='*',color='C2',ecolor='C2',ms=0,alpha=1,lw=0,elinewidth=1,zorder=998)
    
    #mapa1_proflat.plot(df_PDF.prof,df_PDF.lat,'.',ms=1,alpha=0.5)
    mapa1_proflat.plot(hypo[2],hypo[1],'*',mfc='C3',markersize=10,mec='k',zorder=999)
    mapa1_proflat.errorbar(hypo[2],hypo[1],yerr=erlatH,xerr=hypo[4],
                   marker='*',color='C3',ecolor='C3',ms=0,alpha=1,lw=0,elinewidth=1,zorder=998)
    mapa1_proflat.plot(NLL[2],NLL[1],marker='*',mfc='C2',markersize=10,mec='k',zorder=999)
    mapa1_proflat.errorbar(NLL[2],NLL[1],yerr=erlatN,xerr=erproN,
                   marker='*',color='C2',ecolor='C2',ms=0,alpha=1,lw=0,elinewidth=1,zorder=998)
    #mapa1_proflon.plot(df_PDF.lon,df_PDF.prof,'.',ms=1,alpha=0.5)
    mapa1_proflon.plot(hypo[0],hypo[2],'*',mfc='C3',markersize=10,mec='k',zorder=999)
    mapa1_proflon.errorbar(hypo[0],hypo[2],yerr=hypo[4],xerr=erlonH,
                   marker='*',color='C3',ecolor='C3',ms=0,alpha=1,lw=0,elinewidth=1,zorder=998)
    mapa1_proflon.plot(NLL[0],NLL[2],marker='*',mfc='C2',markersize=10,mec='k',zorder=999)
    mapa1_proflon.errorbar(NLL[0],NLL[2],yerr=erproN,xerr=erlonN,
                   marker='*',color='C2',ecolor='C2',ms=0,alpha=1,lw=0,elinewidth=1,zorder=998)
    

    mapa1_proflat.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    mapa1_proflon.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    
    data['x'] = data['x'].round(4)
    data['y'] = data['y'].round(4)
    loncrop = data['x'].where((data['x']<mapa1.get_xlim()[1]) & (data['x']>mapa1.get_xlim()[0]),drop=True)
    latcrop = data['y'].where((data['y']<mapa1.get_ylim()[1]) & (data['y']>mapa1.get_ylim()[0]),drop=True)
    
    zcrop_lon = data['z'].where((data['y']==float(np.round(float(latcenter),2))) & (data['x']<mapa1.get_xlim()[1]) & (data['x']>mapa1.get_xlim()[0]),drop=True)
    zcrop_lat = data['z'].where((data['x']==float(np.round(float(loncenter),2))) & (data['y']<mapa1.get_ylim()[1]) & (data['y']>mapa1.get_ylim()[0]),drop=True)
    mapa1_proflon.plot(loncrop,zcrop_lon[0]*-1/1000,c='C1')
    mapa1_proflat.plot(np.array(zcrop_lat).squeeze()*-1/1000,latcrop,c='C1')
    
    mapa1.vlines(loncenter,mapa1.get_ylim()[0],mapa1.get_ylim()[1],colors='C1',alpha=0.5)
    mapa1.hlines(latcenter,mapa1.get_xlim()[0],mapa1.get_xlim()[1],colors='C1',alpha=0.5)
     
    volnet = gdb.get_metadata_wws(volcan='*')
    volnet = volnet[volnet.tipo=='SISMOLOGICA']
    estas,fases =[],[]
    for i in range(0,len(cat[0].picks)):
        esta = cat[0].picks[i].waveform_id.station_code
        fase = cat[0].picks[i].phase_hint
        estas.append(esta[:-1])
        fases.append(fase)
    estadf = pd.DataFrame([estas,fases]).T
    estadf = estadf[~estadf[0].isin(blacklist)]

    estaP = estadf[estadf[1]=='P']
    estaP = volnet[volnet.codcorto.isin(estaP[0])]
    for index,row in estaP.iterrows():
        mapa1.plot(row.longitud,row.latitud,'o',ms=6,color='C3',fillstyle='top',label='Fase P')
    #estaciones con fase S
    estaS = estadf[estadf[1]=='S']
    estaS = volnet[volnet.codcorto.isin(estaS[0])]
    for index,row in estaS.iterrows():
        mapa1.plot(row.longitud,row.latitud,'o',ms=6,color='C3',fillstyle='full',label='Fases P y S')
    
    mapa1_proflon.yaxis.set_label_position("right")
    mapa1_proflon.yaxis.tick_right()
    
    mapa1_proflat.xaxis.set_label_position("top")
    mapa1_proflat.xaxis.tick_top()
    
    handles, labels = mapa1.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    leg = mapa1.legend(by_label.values(), by_label.keys(),frameon=True)
    i=0
    for lh in leg.legendHandles: 
        if i==0:
            lh._legmarker.set_markersize(10)
        lh._legmarker.set_alpha(1)
        i+=1
        
    mapa1_proflat.set_ylabel('Latitud (°)')
    mapa1_proflat.set_xlabel('Profundidad (km)')
    mapa1_proflon.set_ylabel('Profundidad (km)')
    mapa1_proflon.set_xlabel('Longitud (°)')
    
    mapa1_texto = fig.add_subplot(mapaloc[2])
    mapa1_texto.axis('off')
    mapa1_texto.table(cellText=[['Volcán',voldf.nombre.iloc[0]],['Fecha ',str(evsel.fecha.iloc[0])[0:19]],
                                        ['ML',evsel.ml.iloc[0]],
                                        ['Fases P',len(estaP)],['Fases S',len(estaS)]],colWidths=[0.3,0.75],
                              bbox = [0,0,1,1],zorder=1, rasterized=True) 
  
    trazaax = fig.add_subplot(grilla_externa[1]) 
    
    
    fases = gdb.get_fasesLoc(evsel.idevento.iloc[0], evsel.fecha.iloc[0].year)
    fases = fases[~fases.cod.isin(blacklist2)]

    import geopy.distance
    import numpy as np
    
    
    
    
    dists=[]
    for index,row in fases.iterrows():
        
        estadf = gdb.get_metadata_wws(volcan='*',estadoInst='*',estadoEst='*')
     
        stacod=row.cod
        try:
            stalat = float(estadf[estadf.codcorto==stacod[1:]].latitud.iloc[0])
            stalon = float(estadf[estadf.codcorto==stacod[1:]].longitud.iloc[0])
            stahei = estadf[estadf.codcorto==stacod[1:]].altitud.iloc[0]/1000
         
            evlat = evsel.latitud.iloc[0]
            evlon = evsel.longitud.iloc[0]
            evcoor = (evlat,evlon)
            stacoor = (stalat,stalon)
            
    
            distepi=float(dDICT[stacod[1:]+'Z'])
            print(distepi)
            distprof=stahei+evsel.profundidad_abs.iloc[0]
            distatotal=np.sqrt(distepi**2+distprof**2)
            dists.append(distatotal)
        except:
            dists.append(np.nan)
    fases['dists']=dists
    fases = fases.dropna(subset=['dists'])
    traces =[]
    df_toplot=[]
    import ovdas_SeismicProc_lib as sp
    import ovdas_WWS_lib as WWS
    from datetime import datetime, timedelta
    
    for index,row in fases.iterrows():
        traza = WWS.extraer_signal(estacion=row.cod[1:],inicio=evsel.fecha.iloc[0]-timedelta(seconds=10),fin=evsel.fecha.iloc[0]
                               +timedelta(seconds=evsel.duracion.iloc[0])-timedelta(seconds=1)
                               ,formato_sal='Stream',componente='Z')
        traza = sp.filtrar_traza(traza=traza,tipo="butter",orden=4,fi=1,ff=15)
        traza = traza.merge(fill_value='interpolate')

        
        amps=traza[0].data
        times=traza[0].times('matplotlib')
     
        amps=np.array(amps)
        #amps=amps[amps == amps.astype(float)]
        stap=row.tiempop
        stas=row.tiempos
        
  
       
        amps=amps/abs(max(amps,key=abs))+row.dists
        df_toplot.append([row.cod[1:],row.dists,amps,times,stap,stas,row.dists])
        
    import matplotlib.dates as mdates    
    df_toplot = pd.DataFrame(df_toplot)
    df = df_toplot.sort_values(by=1).reset_index(drop=True)
    pdatas=[]
    ptimes=[]
    for index,row in df.iterrows():
        trazaax.plot(row[3],row[2],linewidth=0.75,color='k',alpha=0.75)
        trazaax.vlines(mdates.epoch2num(row[4]),-0.5+row[6],0.5+row[6],color='C1',lw=2,zorder=999,label='P')   
        trazaax.vlines(mdates.epoch2num(row[5]),-1+row[6],1+row[6],color='C3',lw=2,zorder=999,label='S')
        props = dict(boxstyle='round', facecolor='wheat', alpha=1)
        trazaax.text(max(times),row[6],row[0],ha='center', bbox=props)

    def is_dst(dt=None, timezone="UTC"):
        import pytz
        from datetime import datetime
        if dt is None:
            dt = datetime.utcnow()
        timezone = pytz.timezone(timezone)
        timezone_aware_date = timezone.localize(dt, is_dst=None)
        return timezone_aware_date.tzinfo._dst.seconds != 0
    
    hv=is_dst(timezone='Chile/Continental')
    if hv==True:hvf=-3
    else:hvf=-4
    from datetime import timedelta
    ot=ot-timedelta(hours=hvf)

        
    maxx = 1+(row[6])
    
    trazaax.set_xlim(ot,max(times))


    
    import matplotlib.dates as md
    yfmt = md.DateFormatter('%H:%M:%S')
    trazaax.xaxis.set_major_formatter(yfmt)    
    trazaax.set_ylim(0,maxx)
    trazaax.set_xlabel('Hora (HH:MM:SS)')
    trazaax.set_ylabel('Distancia al origen NLL (km)')
    handles, labels = trazaax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    trazaax.legend(by_label.values(), by_label.keys(),loc=(1.05,0.9))

    
    fig.suptitle('Sistema automático OVV - Sismo '+evsel.tipoevento.iloc[0]+' localizado - '+voldf.nombre.iloc[0]+' - '+str(evsel.fecha.iloc[0])[0:19]+'\n'+
                 'Epicentro (lon,lat): '+str(np.round(hypo[0],2))+'°,'+str(np.round(hypo[1],2))+'° (Hypo71) / '+
                 str(np.round(NLL[0],2))+'°,'+str(np.round(NLL[1],2))+'° (NLL)',fontsize=20)
    


    from PIL import Image
    
    if exec_point == 'local':
        raiz='//172.16.40.102/Monitoreo/ovv2023/'
    elif exec_point == 'server':
        raiz='/mnt/puntocientodos/ovv2023/' 
        
        
    import os
    
    plt.savefig(raiz+'loc_'+str(idev)+'.png')
    jpg = Image.open(raiz+'loc_'+str(idev)+'.png')
    jpg = jpg.convert('RGB')
    jpg.save(raiz+'figs/loc_'+str(idev)+'.jpg') 
    jpg.save(raiz+'last/lastloc.png')     
    
             
    os.remove(raiz+'loc_'+str(idev)+'.png')

    #guardar en subcarpeta figs
    jpg.save('figs/loc_'+str(idev)+'.jpg')
    jpg.save('last/lastloc.png')    
       
    plt.close('all')