# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 09:27:34 2022

@author: sergio.morales
"""

from main import mainnll
import sched, time
s = sched.scheduler(time.time, time.sleep)


def do_something(sc): 
    from datetime import datetime
    fecha=datetime.now()
    print(str(fecha)+ " Monitoreando eventos nuevos clasificados...")
    try:
        mainnll(exec_point='local',debug=False)
    except:
        print('algo pasó :(')        
    s.enter(5, 1, do_something, (sc,))
    
    
    
s.enter(1, 1, do_something, (s,))
s.run()   