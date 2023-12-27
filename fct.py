# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 11:11:48 2023

@author: Gilles Boonen
"""
import numpy as np
import matplotlib.pylab as plt 
from scipy.special import lambertw 

def GetProfile(filename, Shouse=1.0):
    dpm = [0,31,28,31,30,31,30,31,31,30,31,30]
    ConsoData = open(filename, 'r')
    Buf = ConsoData.readlines()
    n_lines = len(Buf)
    dt = 1.
    n_step = int(24. / dt)
    n_days = 365
    E_house = np.zeros((n_days,n_step), dtype=float)
    
    date = np.arange(0,n_days)
    time = np.arange(0,24,dt)

    j = 0 # index for the time
    k = 0 # index for the day
    for i in range(n_lines-1):
        (Buf1, Buf2) = Buf[i+1].split(';')
        (tmp1, tmp2) = Buf1.split(' ')        
        (dd, mm, yy) = tmp1.split('-')
        (h, m) = tmp2.split(':')
        t_float = round(float(h)+float(m)/60.,0)
        j, = np.where(time==t_float)
        k = int(dd)-1+np.sum(dpm[0:int(mm)])
        E_house[k,j] = float(Buf2)
    return time, date, E_house*Shouse

def GetPVProduction(filename, Pnom, fulloutput=False):
    dpm = [31,28,31,30,31,30,31,31,30,31,30]
    PVGISData = open(filename, 'r')
    Buf = PVGISData.readlines()
    n_lines = len(Buf)-11-11
    dt = 1.
    n_step = int(24. / dt)
    n_days = int(n_lines / n_step)
    P_pv = np.zeros((n_days,n_step), dtype=float)
    T2m = np.zeros((n_days,n_step), dtype=float)
    WS10m = np.zeros((n_days,n_step), dtype=float)
    date = np.arange(0,n_days)
    time = np.arange(0,24,dt)
    j = 0 # index for the time
    k = 0 # index for the day
    Ppeak = float(1/Pnom)

    for i in range(n_lines-11):
        (Buf1, Buf2, Buf3, Buf4, Buf5, Buf6, Buf7) = Buf[i+11].split(';')
        (tmp1, tmp2) = Buf1.split(':')
        yy = tmp1[0:4]
        mm = tmp1[4:6]
        dd = tmp1[6:8]
        h = tmp2[0:2]
        m = tmp2[2:4]
        t_float = float(h)
        j, = np.where(time==t_float)
        k = int(dd)-1+int(np.sum(dpm[0:int(mm)-1]))
        P_pv[k,j[0]] = float(Buf2)/1000.
        T2m[k,j[0]] = float(Buf5)
        WS10m[k,j[0]] = float(Buf6)

    if fulloutput: return time, date, P_pv/Ppeak, T2m, WS10m
    else: return time, date, P_pv/Ppeak

def LoadDurationCurve(t, d, P, color='blue', name='Power', steps=250):
    m = len(d) # the number of lines in P
    n = len(t) # the number of columns in P
    dt = 24. / n
    
    Psteps = np.linspace(0., np.amax(P)*1.1, num=steps) # the power considered in the cumulative sum
    Tau = np.zeros(steps, dtype=float) # the number of operating hours
    Pflat = np.reshape(P, m*n)
    for i in range(steps, 0, -1): Tau[i-1] = (np.ma.masked_greater_equal(Pflat, Psteps[i-1])).mask.sum()*dt
    E_tot = 0
    for i in range(m): E_tot += np.trapz(P[i,:], t)
        
    P_max = np.amax(Pflat)
    tau_e = E_tot / P_max
    load_factor = tau_e / (m*24.)
    P_avg = P_max * load_factor

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.close(fig)
    ax.fill_between(Tau, Psteps, 0., facecolor=color, alpha=.5, interpolate=True, label=u'%s= %.0f MWh/y' %(name,E_tot/1000.))
    ax.set_xlim(0, m*24.)
    ax.set_ylim(bottom=0)
    ax.set_ylabel('Power [kW]')
    ax.set_xlabel('Time [h/year]')
    ax.grid(True)
    ax.annotate(r'$P_{peak}$ = %.2f kW ($P_{avg}$ = %.2f kW)' %(P_max, P_avg), xy=(0, P_max), xytext=(tau_e/3, P_max),
            arrowprops=dict(facecolor='black', shrink=0.05))
    ax.text(tau_e, P_max/3*2, r'$\tau_e$ = %.0f h/y (load factor= %.2f)' %(tau_e, load_factor))
    ax.legend()
    plt.show()
    return fig

def CalculFacture(P_el_total, tau_pv, buyCost, sellCost):
    d = 365
    h = 24
    Eres = 0    #Élec tirée du réseau
    Einj = 0    #Élec injectée sur le réseau
    Eauto = 0   #Élec autoconsommée
    for i in range(d):
        for j in range(h):
            dif = P_el_total[i,j] - tau_pv[i,j]
            if dif > 0: #E consommée > E produite
                Eres += dif
                Eauto += tau_pv[i,j]
            elif dif <= 0: #E produite > E consommée
                Eauto += P_el_total[i,j]
                Einj += -dif
    print('\nElec achetée : {}'.format(int(Eres)),' kWh\nElec revendue : {}'.format(int(Einj)),' kWh\nElec autoconsommée : {}'.format(int(Eauto)),' kWh\n')
    cost = Eres/1000 * buyCost - Einj/1000 * sellCost
    print('Facture : {}'.format(int(cost)),' €')
    return cost
  
def Transmission (Theta1):
    n1 = 1          #indice réfraction air
    n2 = 1.5        #indice réfraction verre 
    e = 3           #epaisseur verre
    ke = 5 # Je trouve pas    
    Theta2= np.arcsin((n1*np.sin(Theta1))/n2) #Formule réfraction)
    rho = 0.5*((np.sin(Theta2 - Theta1)**2/np.sin(Theta2 + Theta1)**2)+(np.tan(Theta2 - Theta1)**2/np.tan(Theta2 + Theta1)**2))
    Tau_r = (1 - rho)/(1 + rho) #réflexion
    Tau_a = np.exp(-ke*e/np.cos(Theta2))    #absorption
    Tau = Tau_r * Tau_a         #Trnasmission complète

    return Tau

def GlassAbsorbtion(S, D, theta, n1=1., n2=1.5, e=3., k=3e-3, tau_t=0.95):
    theta_2 = np.degrees(np.arcsin(np.sin(np.radians(theta))*n1/n2))
    rho = np.where(theta != 0, .5*((np.sin(np.radians(theta_2)-np.radians(theta)))**2/(np.sin(np.radians(theta_2)+np.radians(theta)))**2 +(np.tan(np.radians(theta)-np.radians(theta)))**2/(np.tan(np.radians(theta_2)+np.radians(theta)))**2), 0.)
    tau_r = (1.-rho)/(1.+rho)
    tau_a = np.exp(-k*e/np.cos(np.radians(theta_2)))
    tau = np.where(theta<90, tau_a*tau_r, 0.)
    return S*tau+D*tau_t

def PV_cell(V, Gabs, Tamb, NOCT, Rs , Rsh, Id0, T0, Iphi0, G0, n, Eg):
    q = 1.602e-19
    k = 1.38e-23
    Tcell = 273.15 + Tamb + Gabs*(NOCT - 20.)/800.
    Vt = k*Tcell/q
    Id = Id0*(T0/Tcell)**3 * np.exp(q*Eg/n/k*(1./T0 - 1./Tcell))
    Iphi = Iphi0*Gabs/G0
    temp = lambertw(Id*Rs/(n*Vt*(1.+Rs/Rsh))*np.exp(V/n/Vt*(1.-Rs/(Rs+Rsh))+(Iphi+Id)*Rs/(n*Vt*(1.+Rs/Rsh))))
    I = (Iphi + Id - V/Rsh)/(1.+Rs/Rsh) - n*Vt/Rs*temp.real
    return -I*V

def PV_panel(Gabs, Tamb, Scell=15.6*15.6, NOCT=46, Rs = 1., Rsh = 10000., Id0 = 1e-12, T0 = 25+273.15, Iphi0 = 35e-3, G0=1000., n=1., Eg=1.12, ncs=70, ncp=1):
    m = len(Gabs)
    Pm=np.zeros(m, dtype=float)
    Vm=np.zeros(m, dtype=float)
    for i in range(m):
        res = minimize_scalar(PV_cell, args=(Gabs[i], Tamb[i], NOCT, Rs, Rsh, Id0, T0, Iphi0, G0, n, Eg), method='brent', tol=None)
        Pm[i] = -Scell*PV_cell(res.x, Gabs[i], Tamb[i], NOCT, Rs, Rsh, Id0, T0, Iphi0, G0, n, Eg)*ncs*ncp
        Vm[i] = res.x * ncs
    return Pm, Vm

