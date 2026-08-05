
"""
Adiab_CEE 

"""

import numpy as np
import os

logMenv_list=np.linspace(-4,-1,31)
eta_list=np.linspace(0,0.1,11)

logMi_list_HG= np.linspace(np.log10(1.2),np.log10(110),100)
logRi_list_HG= np.linspace(np.log10(2),np.log10(2500),200)

logMi_list_RGB=np.linspace(np.log10(0.2),np.log10(21),100)
logRi_list_RGB=np.linspace(np.log10(0.9),np.log10(1200),200)

logMi_list_AGB=np.linspace(np.log10(0.5),np.log10(21),100)
logRi_list_AGB=np.linspace(np.log10(15),np.log10(1300),200)

logMi_list_massive=np.linspace(np.log10(19),np.log10(110),100)
logRi_list_massive=np.linspace(np.log10(15),np.log10(3300),200)

def RecCEE_solver(Mi,Ri,Menv,stage,core,use_grid):

    logMi=np.log10(Mi);logRi=np.log10(Ri);logMenv=np.log10(Menv)

    if Menv==0.1:
        logMenv=np.log10(0.1-1e-10)
    
    if stage=='HG':
        logMi_list=logMi_list_HG
        logRi_list=logRi_list_HG
    if stage=='RGB':
        logMi_list=logMi_list_RGB
        logRi_list=logRi_list_RGB
    if stage=='AGB':
        logMi_list=logMi_list_AGB
        logRi_list=logRi_list_AGB
    if stage=='massive':
        logMi_list=logMi_list_massive
        logRi_list=logRi_list_massive
    
    #----------determine x-------------#
    for x in range(100):
        if logMi<logMi_list[x]:
            logMi_left =logMi_list[x-1]
            logMi_right=logMi_list[x]
            break
    
    for y in range(200):
        if logRi<logRi_list[y]:
            logRi_left =logRi_list[y-1]
            logRi_right=logRi_list[y]
            break
    
    for z in range(31):
        if logMenv<logMenv_list[z]:
            logMenv_left =logMenv_list[z-1]
            logMenv_right=logMenv_list[z]
            break
    
    xd=(logMi-logMi_left)/(logMi_right-logMi_left)
    yd=(logRi-logRi_left)/(logRi_right-logRi_left)
    zd=(logMenv-logMenv_left)/(logMenv_right-logMenv_left)
    
    ##------------interpolation--------------##
    if core=='core1':
        loglam0_name='loglam0_grid_'+stage+'_0.01.npy'
        F1_name='F1_grid_'+stage+'_0.01_'+use_grid+'.npy'
        F2_name='F2_grid_'+stage+'_0.01_'+use_grid+'.npy'
    if core=='core2':
        loglam0_name='loglam0_grid_'+stage+'_0.1.npy'
        F1_name='F1_grid_'+stage+'_0.1_'+use_grid+'.npy'
        F2_name='F2_grid_'+stage+'_0.1_'+use_grid+'.npy'
        
    loglam0_list=np.load('data/'+loglam0_name)
    F1_list=np.load('data/'+F1_name)
    F2_list=np.load('data/'+F2_name)
        
    #----------on lam0----------#
    lam_00=loglam0_list[x-1][y-1]
    lam_01=loglam0_list[x-1][y]
    lam_10=loglam0_list[x][y-1]
    lam_11=loglam0_list[x][y]
    
    u_0=lam_00*(1-xd)+lam_10*xd
    u_1=lam_01*(1-xd)+lam_11*xd
    loglam0_3d=u_0*(1-yd)+u_1*yd
    lam0_3d=10**loglam0_3d
    #print(loglam0_3d)
    
    #----------on F1----------#
    
    F1_000=F1_list[x-1][y-1][z-1]
    F1_001=F1_list[x-1][y-1][z]
    F1_010=F1_list[x-1][y][z-1]
    F1_011=F1_list[x-1][y][z]
    F1_100=F1_list[x][y-1][z-1]
    F1_101=F1_list[x][y-1][z]
    F1_110=F1_list[x][y][z-1]
    F1_111=F1_list[x][y][z]
    #print(F1_000,F1_111)
    
    v_00=F1_000*(1-xd)+F1_100*xd
    v_01=F1_001*(1-xd)+F1_101*xd
    v_10=F1_010*(1-xd)+F1_110*xd
    v_11=F1_011*(1-xd)+F1_111*xd
    
    v_0=v_00*(1-yd)+v_10*yd
    v_1=v_01*(1-yd)+v_11*yd
    F1_3d=v_0*(1-zd)+v_1*zd
    
    if F1_3d<1:
        F1_3d=1.0
    
    #print(F1_3d)
    
    #----------on F2----------#
    
    F2_000=F2_list[x-1][y-1][z-1]
    F2_001=F2_list[x-1][y-1][z]
    F2_010=F2_list[x-1][y][z-1]
    F2_011=F2_list[x-1][y][z]
    F2_100=F2_list[x][y-1][z-1]
    F2_101=F2_list[x][y-1][z]
    F2_110=F2_list[x][y][z-1]
    F2_111=F2_list[x][y][z]
    #print(F2_000,F2_111)
    
    w_00=F2_000*(1-xd)+F2_100*xd
    w_01=F2_001*(1-xd)+F2_101*xd
    w_10=F2_010*(1-xd)+F2_110*xd
    w_11=F2_011*(1-xd)+F2_111*xd
    
    w_0=w_00*(1-yd)+w_10*yd
    w_1=w_01*(1-yd)+w_11*yd
    F2_3d=w_0*(1-zd)+w_1*zd
    
    if F2_3d<1:
        F2_3d=1.0
    
    #print(F2_3d)

    return [lam0_3d,F1_3d,F2_3d]

def Rfad_solver(Mi,Ri,Menv,stage,core):

    logMi=np.log10(Mi);logRi=np.log10(Ri);eta=Menv/Mi

    if eta==0.1:
        eta=0.1-1e-10
    
    if stage=='HG':
        logMi_list=logMi_list_HG
        logRi_list=logRi_list_HG
    if stage=='RGB':
        logMi_list=logMi_list_RGB
        logRi_list=logRi_list_RGB
    if stage=='AGB':
        logMi_list=logMi_list_AGB
        logRi_list=logRi_list_AGB
    if stage=='massive':
        logMi_list=logMi_list_massive
        logRi_list=logRi_list_massive
    
    #----------determine x-------------#
    for x in range(100):
        if logMi<logMi_list[x]:
            logMi_left =logMi_list[x-1]
            logMi_right=logMi_list[x]
            break
    
    for y in range(200):
        if logRi<logRi_list[y]:
            logRi_left =logRi_list[y-1]
            logRi_right=logRi_list[y]
            break
    
    for z in range(11):
        if eta<eta_list[z]:
            eta_left =eta_list[z-1]
            eta_right=eta_list[z]
            break
    
    xd=(logMi-logMi_left)/(logMi_right-logMi_left)
    yd=(logRi-logRi_left)/(logRi_right-logRi_left)
    zd=(eta-eta_left)/(eta_right-eta_left)
    
    #print(xd,yd,zd)
    
    ##------------interpolation--------------##
    if core=='core1':
        Rf_name='Rf_grid_'+stage+'_0.01'+'.npy'
    if core=='core2':
        Rf_name='Rf_grid_'+stage+'_0.1'+'.npy'

    Rf_list=np.load('data/'+Rf_name)
    
    #----------on Rf-ad----------#
    
    Rf_000=Rf_list[x-1][y-1][z-1]
    Rf_001=Rf_list[x-1][y-1][z]
    Rf_010=Rf_list[x-1][y][z-1]
    Rf_011=Rf_list[x-1][y][z]
    Rf_100=Rf_list[x][y-1][z-1]
    Rf_101=Rf_list[x][y-1][z]
    Rf_110=Rf_list[x][y][z-1]
    Rf_111=Rf_list[x][y][z]
    #print(F1_000,F1_111)
    
    v_00=Rf_000*(1-xd)+Rf_100*xd
    v_01=Rf_001*(1-xd)+Rf_101*xd
    v_10=Rf_010*(1-xd)+Rf_110*xd
    v_11=Rf_011*(1-xd)+Rf_111*xd
    
    v_0=v_00*(1-yd)+v_10*yd
    v_1=v_01*(1-yd)+v_11*yd
    Rf_3d=v_0*(1-zd)+v_1*zd
    
    #print(F1_3d)

    return Rf_3d

# module info
__version__ = "1.0.0"
__author__ = "Lifu Zhang"

# function out put
__all__ = ['RecCEE_solver', 'Rfad_solver']
