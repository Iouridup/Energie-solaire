from scipy.special import lambertw
import numpy as np
import matplotlib.pylab as plt
from ipywidgets import interact, FloatSlider, fixed


Rs = 1.
Rsh = 1000.
T0 = 25. + 273.15
Ta = 35 + 273.15
NOCT = 48.  #Normal operation cell T
Id0 = 10e-13
Iphi0 = 35e-3 #we suppose a radiation of 1000 W/m²
k = 1.38e-23
q = 1.602e-19
Tc = 298.15
Vt = k*Tc/q
Eg = 1.12
n=1.
G=1000.

V = np.linspace(0,.7, num=100)
I = np.zeros(100, dtype=float)

P = np.zeros(100, dtype=float)
for i in range(100):
    Tc = Ta + G*(NOCT-20.)/800.
    Iphi = Iphi0*G/1000.
    Id = Id0*(T0/Tc)**3 * np.exp(q*Eg/n/k*(1./T0 - 1./Tc))
    temp = lambertw(Id*Rs/(n*Vt*(1.+Rs/Rsh))*np.exp(V[i]/n/Vt*(1.-Rs/(Rs+Rsh))+(Iphi+Id)*Rs/(n*Vt*(1.+Rs/Rsh))), 0)
    I[i] = (Iphi + Id - V[i]/Rsh)/(1.+Rs/Rsh) - n*Vt/Rs*temp.real
    P[i] = I[i]*V[i]

Pm = np.amax(P)
Vm = V[np.where(P==Pm)[0][0]]
temp = lambertw(Id*Rs/(n*Vt*(1.+Rs/Rsh))*np.exp(Vm/n/Vt*(1.-Rs/(Rs+Rsh))+(Iphi+Id)*Rs/(n*Vt*(1.+Rs/Rsh))))
Im = (Iphi + Id - Vm/Rsh)/(1.+Rs/Rsh) - n*Vt/Rs*temp.real
fig, ax = plt.subplots(figsize=(12,8))
ax.plot(V,I,'k-', label='Current')
ax.plot(V,P,'r-', label='Power (W)')
ax.axvline(x=Vm, color='k', linestyle=':')
ax.axhline(y=Im, color='k', linestyle=':')
ax.set_xlim(0,0.7)
ax.set_ylim(0,40e-3)
ax.set_xlabel('Voltage (V)')
ax.set_ylabel('Current (A/cm2)')
ax.grid(True)
#ax.legend()

epsilon=Pm*12*12/(G*12*12/10000)

print('A photocell generates %.2f A and %.2f W' %(Im*12*12, Pm*12*12))
print('the efficiency is', epsilon*100, '%')
print('Un PV de 72 cellules produira %.2f A and %.2f W' %(Im*12*12*72, Pm*12*12*72))


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

a = PV_cell
print("PV cell vaut : "a)