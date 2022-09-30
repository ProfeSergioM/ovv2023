# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 12:56:46 2022

@author: sergio.morales
"""
def mainnll(exec_point,debug):
    from lib.syslib import addPath
    from lib.dblib import getLastEv,saveNLLcsv
    from lib.nlllib import doNLL
    from lib.plotlib import plotNLLinfo
    from lib.maillib import sendMail
    
    addPath(exec_point)
    reloc, ev= getLastEv(debug) #obtiene el último evento localizado
    if reloc==True:
        NLL,estadf,hypo,voldf = doNLL(ev,exec_point,fini='',ffin='')
        saveNLLcsv(voldf, ev, NLL, exec_point)
        plotNLLinfo(ev,exec_point)
        sendMail(voldf, ev)