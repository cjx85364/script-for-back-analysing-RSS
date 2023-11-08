#!/usr/bin/env python
# coding: utf-8
# %%

# this file is used to test the best cov



import numpy as np

import copy
# create random cohesion and friction angle couples 
def sample(c,fai,cov,size):
    cohesion=np.random.normal(c,c*cov,size=size)#cohesion
    fai=np.random.normal(fai,fai*0.15,size=size)#fai
    for i in range(len(cohesion)):
        while cohesion[i] <0:
            cohesion[i]=np.random.normal(c,c*cov)
    for i in range(len(fai)):    
        while fai[i] <0:
            fai[i]=np.random.normal(fai,fai*0.15)
    fai=np.reshape(fai,(len(fai),1))
    cohesion=np.reshape(cohesion,(len(cohesion),1))
    return cohesion, fai
# runing a directly Newmark solution  
def calculate_pro(c,phi,thi,slope,den,PGA):
    tan_phi=np.tan(phi*np.pi/180.)
    FS = c/(np.sin(slope*np.pi/180.)*den*thi)+tan_phi/np.tan(slope*np.pi/180.)
    ky=(FS-1)*np.sin(slope*np.pi/180.)#unit/g 
    cuntr=sum(1 for x in ky if x > PGA)
    cuntr = 1 if cuntr < 1 else cuntr
    mean_phi=np.sum(phi[ky>PGA])/cuntr
    mean_C=np.sum(c[ky>PGA])/cuntr
    pro_FS = mean_C/(np.sin(slope*np.pi/180.)*den*thi)+np.tan(mean_phi*np.pi/180.)/np.tan(slope*np.pi/180.)
    pro_ky=(pro_FS-1)*np.sin(slope*np.pi/180.)
    if pro_ky>PGA: 
        soil_probability_of_failure=0
    else:
        soil_probability_of_failure=1
    return soil_probability_of_failure
# give the staring point according to lithology 
def calculate_all(cl,thi,slope,pga,cov):
    pro=np.random.rand(len(pga))
    c_granit,phi_granit=sample(40,44,cov,5000)#26.5
    c_phy,phi_phy=sample(27,35,cov,5000)#27.8
    c_diorite,phi_diorite=sample(40,43,cov,5000)#27
    c_and,phi_and=sample(39,42,cov,5000)#26.2
    c_anor,phi_anor=sample(40,42,cov,5000)#27
    c_basalt,phi_basalt=sample(45,45,cov,5000)#29
    c_lime,phi_lime=sample(34,27.8,cov,5000)#40
    c_sand,phi_sand=sample(29,39,cov,5000)#26
    for i in range(len(cl)):
        if cl[i]=="granite":
            pro[i]=calculate_pro(c_granit,phi_granit,thi,slope[i],26.5,pga[i])
        elif cl[i]=="phyllite":
            pro[i]=calculate_pro(c_phy,phi_phy,thi,slope[i],27.8,pga[i])
        elif cl[i]=="diorite":
            pro[i]=calculate_pro(c_diorite,phi_diorite,thi,slope[i],27,pga[i])
        elif cl[i]=="and":
            pro[i]=calculate_pro(c_and,phi_and,thi,slope[i],26.2,pga[i])
        elif cl[i]=="anor":
            pro[i]=calculate_pro(c_anor,phi_anor,thi,slope[i],27,pga[i])
        elif cl[i]=="basalt":
            pro[i]=calculate_pro(c_basalt,phi_basalt,thi,slope[i],29,pga[i])
        elif cl[i]=="lime":
            pro[i]=calculate_pro(c_lime,phi_lime,thi,slope[i],27.8,pga[i])     
        elif cl[i]=="sand":
            pro[i]=calculate_pro(c_sand,phi_sand,thi,slope[i],26,pga[i])
        else:
            pro[i]=0
    return pro
#calculate the TNR by using the a set of cov values
def cov_cal(start,end,interal,cl,t,slope,pga):
    point = np.arange(start, end, interal)
    a=np.random.rand(len(point))
    for i in range(len(point)):
        pro=calculate_all(cl,t,slope,pga,point[i])
        fn = sum(pro == 1)
        tn = sum(pro == 0)
        a[i]=tn/(fn+tn)
    return a