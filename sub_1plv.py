# -*- coding: utf-8 -*-
"""
@author: C
"""

import CoolProp.CoolProp as CP
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##DATI
T_br_in=150+273.15 #Kelvin
T_br_out_min=70+273.15 #Kelvin
Tml=(T_br_in-T_br_out_min)/np.log((T_br_in/T_br_out_min))
p_br= 0.5e6 #Pascal
m_dot_br=200 #kg/s
br_f='water'

T_k=40+273.15 #Kelvin
eta_pump=0.7
eta_turbine=0.82

tau_pp=5 #Kelvin


## POTENZA TERMICA MASSIMA ESTRAIBILE
h_br_1=CP.PropsSI('H','P',p_br,'T',T_br_in,br_f) #J/kg
h_br_2=CP.PropsSI('H','P',p_br,'T',T_br_out_min,br_f) #J/kg
Q_th_max=m_dot_br*(h_br_1-h_br_2) #Watt

#PROFLIO TERMICO FLUIDO GEOTERMICO
n=100
T_br=np.linspace(T_br_out_min,T_br_in,n)
h_br=CP.PropsSI('H','P',p_br,'T',T_br,br_f)
Q=m_dot_br*(h_br-h_br_2)

#CREAZIONE DATASET FLUIDI
wf_all=CP.get_global_param_string("FluidsList").split(',')
wf_data=pd.DataFrame(index=wf_all, columns=["Pcrit","Tcrit"])
label=wf_data.columns
for wf in wf_all:
   wf_data.loc[wf,wf_data.columns[0]]=CP.PropsSI('Pcrit',wf)
   wf_data.loc[wf,wf_data.columns[1]]=CP.PropsSI('Tcrit',wf)

#SELEZIONE IN BASE ALLA TEMPERATURA CRITICA
iclean=wf_data[label[1]]>T_br_in
wf_data_clean=wf_data[iclean]

def cycle_1plv(p,p_k,fl,T_k,eta_pump,eta_turbine,Q_th_max):
    h1=CP.PropsSI('H','Q',0,'T',T_k,fl)
    s1=CP.PropsSI('S','Q',0,'T',T_k,fl)
    h2is=CP.PropsSI('H','S',s1,'P',p,fl)
    h2=h1+(h2is-h1)/eta_pump
    h3=CP.PropsSI('H','Q',1,'P',p,fl)
    s3=CP.PropsSI('S','Q',1,'P',p,fl)
    h4is=CP.PropsSI('H','S',s3,'P',p_k,fl)
    h4=h3-(h3-h4is)*eta_turbine
    m_dot_fl=Q_th_max/(h3-h2)
    eta=((h3-h4)-(h2-h1))/(h3-h2)
    return m_dot_fl,h1,h2,h3,h4,eta
        

def HEX_profile(fl,T_k,T_br,eta_pump,eta_turbine,Q,Q_th_max):
    p_k=CP.PropsSI('P','T',T_k,'Q',0,fl)
    p_v=CP.PropsSI('Pcrit',fl)
    parr=np.geomspace(p_k*1.005,p_v*0.95,300)
    HEX_data=pd.DataFrame(index=parr, columns=['mdot','eta','pk','delta_min','Torc','h1','h2','h3','h4','h23'])
    for p in parr:
        m_dot,h1,h2,h3,h4,eta=cycle_1plv(p,p_k,fl,T_k,eta_pump,eta_turbine, Q_th_max)
        h_23=h2+Q/m_dot
        T_orc = CP.PropsSI('T','P', p,'H',h_23,fl)
        delta=T_br-T_orc
        delta_min=min(delta)
        HEX_data.loc[p,'Torc']=T_orc
        HEX_data.loc[p,'mdot']=m_dot
        HEX_data.loc[p,'eta']=eta
        HEX_data.loc[p,'delta_min']=delta_min
        HEX_data.loc[p,'h1']=h1
        HEX_data.loc[p,'h2']=h2
        HEX_data.loc[p,'h3']=h3
        HEX_data.loc[p,'h4']=h4
        HEX_data.loc[p,'pk']=p_k
        HEX_data.loc[p,'h23']=h_23
    return HEX_data
    

data_all=pd.DataFrame(index=wf_data_clean.index, columns=['p','pk','beta','vin','vout','Vr','lc','lt','mdot','tau_pp','eta','Power[MW]'])


for wf in wf_data_clean.index:
    
    dati_wf=HEX_profile(wf,T_k,T_br,eta_pump,eta_turbine,Q,Q_th_max)
    
    best=dati_wf['delta_min']>=tau_pp
    dati_wf=dati_wf[best]
    
    best=dati_wf['eta']==max(dati_wf['eta'])
    dati_wf=dati_wf[best]
    
    data_all.loc[wf,'p']=dati_wf.index[0]/1e5 #bar
    data_all.loc[wf,'pk']=dati_wf['pk'].iloc[0]/1e5 #bar
    data_all.loc[wf,'beta']=dati_wf.index[0]/dati_wf['pk'].iloc[0]
    data_all.loc[wf,'vin']=1/(CP.PropsSI('D','H',dati_wf['h3'].iloc[0],'P',dati_wf['pk'].iloc[0],wf))
    data_all.loc[wf,'vout']=1/(CP.PropsSI('D','H',dati_wf['h4'].iloc[0],'P',dati_wf.index[0],wf))
    data_all.loc[wf,'Vr']=data_all.loc[wf,'vin']/data_all.loc[wf,'vout']
    data_all.loc[wf,'lc']=(dati_wf['h2'].iloc[0]-dati_wf['h1'].iloc[0])/1e3 #kJ/kg
    data_all.loc[wf,'lt']=(dati_wf['h3'].iloc[0]-dati_wf['h4'].iloc[0])/1e3 #kJ/kg
    data_all.loc[wf,'mdot']=dati_wf['mdot'].iloc[0] #kg/s
    data_all.loc[wf,'tau_pp']=dati_wf['delta_min'].iloc[0]
    data_all.loc[wf,'eta']=dati_wf['eta'].iloc[0]
    data_all.loc[wf,'Power[MW]']=(data_all.loc[wf,'lt']-data_all.loc[wf,'lc'])*data_all.loc[wf,'mdot']/1e3
    
    # hi=dati_wf['h23'].iloc[0]
    # p=dati_wf.index[0]
    # si=CP.PropsSI('S','H',hi,'P',p,wf)
    # T=CP.PropsSI('T','S',si,'P',p,wf)
    # Tint=np.linspace(273.15+20, wf_data_clean.loc[wf,'Tcrit'])
    # sls=CP.PropsSI('S','Q',0,'T',Tint,wf)
    # svs=CP.PropsSI('S','Q',1,'T',Tint,wf)
    # plt.title('['+wf+']'+' '+'T-s')
    # plt.plot(sls/1e3,Tint-273.15,linestyle=':',color='0.4',linewidth=1.0,label='Saturation Curve')
    # plt.plot(svs/1e3,Tint-273.15,linestyle=':',color='0.4',linewidth=1.0)
    # plt.plot(si/1e3,T-273.15,'green',label=wf)
    # plt.plot(si/1e3,T_br-273.15,'red',label='Brine')
    # plt.xlabel('Entropia specifica [kJ/(kgK)]')
    # plt.ylabel('Temperatura [°C]')
    # plt.ylim(30,170)
    # plt.ylim()
    # plt.legend()
    # plt.grid()
    # # plt.savefig(f'{i}.pdf',format='pdf',bbox_inches='tight')
    # plt.show()
    # plt.close()

    # plt.title('['+wf+']'+' '+'T-Q')
    # plt.plot(Q/Q_th_max*100,T_br-273.15,'red',label='Brine')
    # plt.plot(Q/Q_th_max*100,dati_wf['Torc'].iloc[0]-273.15,'blue',label=wf)
    # plt.xlabel('Frazione di Potenza termica [%]')
    # plt.ylabel('Temperatura [°C]')
    # plt.legend()
    # plt.grid()
    # plt.show()
    # plt.close()



