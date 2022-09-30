# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:18:43 2022

@author: sergio.morales
"""

def addPath(exec_point):
    import locale
    locale.setlocale(locale.LC_ALL, 'en_US.utf8')
    import os, sys
    sys.path.insert(0, os.path.abspath('..'))
    import sys
    if exec_point == 'local':
        sys.path.append('//172.16.40.10/sismologia/pyovdas_lib/')
    elif exec_point == 'server':
        sys.path.append('/mnt/puntodiez/pyovdas_lib/')
    print('Path '+exec_point)