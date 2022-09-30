# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 16:19:43 2022

@author: sergio.morales
"""

import os
import pickle
# Gmail API utils
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
# for encoding/decoding messages in base64
#from base64 import urlsafe_b64decode, urlsafe_b64encode
# for dealing with attachement MIME types
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.image import MIMEImage
from email.mime.audio import MIMEAudio
from email.mime.base import MIMEBase
from mimetypes import guess_type as guess_mime_type
#import textwrap

# Request all access (permission to read/send/receive emails, manage the inbox, and more)
SCOPES = ['https://mail.google.com/']
our_email = 'ovvearlywarning@gmail.com'

def gmail_authenticate():
    creds = None
    # the file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first time
    if os.path.exists("keys/token.pickle"):
        with open("keys/token.pickle", "rb") as token:
            creds = pickle.load(token)
    # if there are no (valid) credentials availablle, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file('keys/gmailapi.json', SCOPES)
            creds = flow.run_local_server(port=0)
        # save the credentials for the next run
        with open("keys/token.pickle", "wb") as token:
            pickle.dump(creds, token)
    return build('gmail', 'v1', credentials=creds)



# Adds the attachment with the given filename to the given message
def add_attachment(message, filename):
    content_type, encoding = guess_mime_type(filename)
    if content_type is None or encoding is not None:
        content_type = 'application/octet-stream'
    main_type, sub_type = content_type.split('/', 1)
    if main_type == 'text':
        fp = open(filename, 'rb')
        msg = MIMEText(fp.read().decode(), _subtype=sub_type)
        fp.close()
    elif main_type == 'image':
        fp = open(filename, 'rb')
        msg = MIMEImage(fp.read(), _subtype=sub_type)
        fp.close()
    elif main_type == 'audio':
        fp = open(filename, 'rb')
        msg = MIMEAudio(fp.read(), _subtype=sub_type)
        fp.close()
    else:
        fp = open(filename, 'rb')
        msg = MIMEBase(main_type, sub_type)
        msg.set_payload(fp.read())
        fp.close()
    filename = os.path.basename(filename)
    msg.add_header('Content-Disposition', 'attachment', filename=filename)
    message.attach(msg)
    
    
def build_message(destination, obj, body, attachments=[]):
    from base64 import urlsafe_b64encode
    cuerpo = MIMEText(body,'plain')
    
    
    if not attachments: # no attachments given
        message = cuerpo
        message['to'] = destination
        message['from'] = our_email
        message['subject'] = obj
    else:
        message = MIMEMultipart()
        message['to'] = destination
        message['from'] = our_email
        message['subject'] = obj
        message.attach(cuerpo)
        for filename in attachments:
            add_attachment(message, filename)
    return {'raw': urlsafe_b64encode(message.as_bytes()).decode()}


def send_message(service, destination, obj, body, attachments=[]):
    return service.users().messages().send(
      userId="me",
      
      body=build_message(destination, obj, body, attachments)
    ).execute()


def sendMail(voldf,evsel,tipomail,debug):
    import numpy as np
    #volcan=voldf.nombre_db.iloc[0]
    idev = evsel.idevento.iloc[0]
    amp= evsel.amplitud_ums.iloc[0]
    frec= evsel.frecuencia.iloc[0]
    tipoev = evsel.tipoevento.iloc[0]
    prof = evsel.profundidad_abs.iloc[0]
    lat = evsel.latitud.iloc[0]
    lon = evsel.longitud.iloc[0]
    ML= evsel.ml.iloc[0]
    volcan_real = voldf.nombre.iloc[0]
    tipo = voldf.vol_tipo.iloc[0]
    if debug==True:
        lista=['sergiomoralesmendez@gmail.com','sergio.morales@sernageomin.cl']
    else:
        lista = ['sergiomoralesmendez@gmail.com',
                        'sergio.morales@sernageomin.cl',
                        'luis.franco@sernageomin.cl',
                        'fernando.gil@sernageomin.cl',
                        'carlos.cardona@sernageomin.cl',
                        'paola.pena@sernageomin.cl',
                        'oscar.valderrama@sernageomin.cl',
                        'jonathan.quijada@sernageomin.cl',
                        'juan.sanmartin@sernageomin.cl',
                        'monitoreo1@sernageomin.cl',
                        'monitoreo2@sernageomin.cl',
                        'monitoreo3@sernageomin.cl'
                        ]

  
    if tipomail=='REAV':
        asunto="PRUEBA - REAV "+tipoev+" ML="+str(ML)+" "+tipo+" "+volcan_real+" - PRUEBA"
        cuerpo = """\nSismo localizado en el {tipov} {volcan}:
        \nTipo Evento: {tipo}
    Magnitud : {ml} (ML)
    Amplitud: {amplitud} um/s (referencia)
    Frecuencia: {frecuencia} Hz (Referencia)
    Profundidad: {profundidad} km bajo nivel del mar
    Latitud: {latitud}°
    Longitud: {longitud}°
    \nSe adjunta Propuesta de REAV.
    \nSaludos.
    \n(Este correo es automático, no responder)
        """
        cuerpo = cuerpo.format(tipo=tipoev,ml=ML,tipov=tipo,volcan=volcan_real,amplitud=amp,frecuencia=frec,profundidad=np.round(prof,1),
                               latitud=np.round(lat,2),longitud=np.round(lon,2))
    elif tipomail=='normal':
        asunto="Sismo tipo "+tipoev+" ML="+str(ML)+" "+tipo+" "+volcan_real
        cuerpo = """\nSismo localizado en el {tipov} {volcan}:
        \nTipo Evento: {tipo}
    Magnitud : {ml} (ML)
    Amplitud: {amplitud} um/s (referencia)
    Frecuencia: {frecuencia} Hz (Referencia)
    Profundidad: {profundidad} km bajo nivel del mar
    Latitud: {latitud}°
    Longitud: {longitud}°
    \nSe adjunta imagen correspondiente a localización con Hypo71 realizada en la sala de monitoreo y localización usando NLL en base a las fases P y S detectadas por los especialistas.
    \nSaludos.
    \n(Este correo es automático, no responder)
        """
        cuerpo = cuerpo.format(tipo=tipoev,ml=ML,tipov=tipo,volcan=volcan_real,amplitud=amp,frecuencia=frec,profundidad=np.round(prof,1),
                               latitud=np.round(lat,2),longitud=np.round(lon,2))        
    for destino in lista:
        print(destino)
        send_message(gmail_authenticate(), destino, asunto,cuerpo,['figs/loc_'+str(idev)+'.jpg']) 
