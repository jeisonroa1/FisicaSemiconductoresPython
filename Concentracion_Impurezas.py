## Concentración de Impurezas
## Jeison Ivan Roa M.

import numpy 	         as np
import scipy.constants   as sc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.figure import Figure

from matplotlib import colors
from Tkinter import *

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

figure = Figure()
canvas = FigureCanvas(figure)

e = sc.e
h = sc.h
k = sc.k 
pi = np.pi 
mo = 9.11e-31
T = np.linspace(100,2000)

"""Los parametros del Silicio, Germanio y el Arsenuro de galio estan organizados de la siguiente forma: Eg(0),m*l,m*t,mll,mlh,Mc; 
para el Arsenuro de Galio la masa efectiva es 0.063*m0"""
GaAs = [1.519,0.063,0.076,0.5,1]
Si = [1.17,0.98,0.19,0.16,0.49,6]
Ge = [0.7437,1.64,0.082,0.04,0.28,8]

"""Parametros del Silicio"""
mdeSi = ((Si[1]*mo)*(Si[2]*mo)**2)**(1.0/3.0)
mdhSi = ((Si[3]*mo)**(3.0/2.0) + (Si[4]*mo)**(3.0/2.0))**(2.0/3.0)

"""Parametros del Germanio"""
mdeGe = ((Ge[1]*mo)*(Ge[2]*mo)**2)**(1.0/3.0)
mdhGe = ((Ge[3]*mo)**(3.0/2.0) + (Ge[4]*mo)**(3.0/2.0))**(2.0/3.0)

"""Parametros del Arsenuro de Galio"""
mdeGaAs = GaAs[1]*mo
mdhGaAs = ((GaAs[2]*mo)**(3.0/2.0) + (GaAs[3]*mo)**(3.0/2.0))**(2.0/3.0)

def niSi(T):
        return np.exp(-Si[0]/(2*(k/e)*T))*(2*Si[5]**(1.0/2.0)*(2*pi*k*T)**(3.0/2.0)*(h)**(-2))*(mdeSi*mdhSi/mo**2)**(3.0/4.0)*10**(-14)

def niGe(T):
        return np.exp(-Ge[0]/(2*(k/e)*T))*(2*Ge[5]**(1.0/2.0)*(2*pi*k*T)**(3.0/2.0)*(h)**(-2))*(mdeGe*mdhGe/mo**2)**(3.0/4.0)*10**(-14)

def niGeAs(T):
	return np.exp(-GaAs[0]/(2*(k/e)*T))*(2*GaAs[4]**(1.0/2.0)*(2*pi*k*T)**(3.0/2.0)*(h)**(-2))*(mdeGaAs*mdhGaAs/mo**2)**(3.0/4.0)*10**(-14)

plt.plot(niSi(T),T, color="b", label=r"$Si$")
plt.plot(niGe(T),T, color="r",  label=r"$Ge$")
plt.plot(niGeAs(T),T, color="g",  label=r"$GaAs$")
plt.legend(loc=4)
plt.xscale('symlog')
plt.title(r"$Temperatura \ \ intrinseca \ \ T_i \ \ como \ \ una \ \ funcion \ \ de \ \ la \ \ concentracion \ \ de \ \ fondo $")
plt.ylabel(r"$Temperatura \ \ intrinseca \ \ T_i \ \ (^{\circ} C)$")
plt.xlabel(r"$Concentracion \ \ de \ \ Impurezas \ \ (cm^{-3})$")
plt.xlim(10**13,2*10**18)
plt.ylim(100,600)
##plt.grid()
plt.show()

