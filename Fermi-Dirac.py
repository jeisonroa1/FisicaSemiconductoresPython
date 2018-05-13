## Fermi-Dirac
## Jeison Ivan Roa Mora

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

fig = plt.figure() 
ax = plt.subplot('111', axisbg='w')


#Constantes utilizadas
k  = sc.k
e = sc.e

#Energia en eV relativa
Eo = 5.0; Emin=0; Emax=10.0;
E = np.linspace(Emin,Emax,1000)

ymin=-0.1; ymax=1.1; xval=0.05; yval =2.3;

line, = ax.plot([], [])

plt.ylabel("$F \ \ (E)$")
plt.xlabel("$Energia \ \ (E) \ \ [e\cdot V ]$")
plt.title('$FUNCION \ \ DE \ \ DISTRIBUCION \ \ FERMI-DIRAC \ \ - \ \ F(E)$')

def fun(Temp, color1, etq):

    #Funcion de distribucion de Fermi-Dirac
    ax.plot(E,((1)/((np.exp((E-Eo)/(k/e*Temp)))+(1.0))),color=plt.cm.gist_rainbow(color1), label=etq)

    
T1 = 0
T2 = 1000
T3 = 2500
T4 = 4500
T5 = 7000
T6 = 10000

col1 = int(200./T6*T1)
col2 = int(200./T6*T2)
col3 = int(200./T6*T3)
col4 = int(200./T6*T4)
col5 = int(200./T6*T5)
col6 = int(200./T6*T6)

fun(T1,col1, r"$F(E)_{T = %d K}$" %T1)
fun(T2,col2, r"$F(E)_{T = %d K}$" %T2)
fun(T3,col3, r"$F(E)_{T = %d K}$" %T3)
fun(T4,col4, r"$F(E)_{T = %d K}$" %T4)
fun(T5,col5, r"$F(E)_{T = %d K}$" %T5)
fun(T6,col6, r"$F(E)_{T = %d K}$" %T6)

legend = plt.legend(loc=1, prop={'size':12.5})
plt.axis([Emin,Emax,ymin,ymax])

def animate(i): 
    y=((1)/((np.exp((E-Eo)/(k/e*(i*50))))+(1.0)))
    line, = ax.plot(E, y, color=plt.cm.gist_rainbow(i))
    return line,

ani = animation.FuncAnimation(fig, animate,  frames=200, interval=50, blit = True, repeat=True)

legend = plt.legend(loc=1)
xano1=((Eo-Emin)/2.0)-0.1
xano2= Eo+((Emax-Eo)/2.0)-0.1
            
ax.annotate(r"$F(E) = \frac{1}{e^\frac{E - E_{f}}{kT}+1}$", xy=(0.5,0.4), color='k',
            fontsize=25,
	    horizontalalignment='left',
            verticalalignment='up',
            )
            
ax.annotate(r"$F(E)_{T = 0K, E<E_{f}} = 1$", xy=(0.5,0.2), color='k',
            fontsize=15,
	    horizontalalignment='left',
            verticalalignment='up',
            )
ax.annotate("$F(E)_{T = 0K, E>E_{f}} = 0$", xy=(0.5,0.1), color='k',
            fontsize=15,
	    horizontalalignment='left',
            verticalalignment='up',
            )
            
ax.annotate("$E_f$", xy=(Eo-0.2,-0.08), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )

ax.annotate("$E > E_f$", xy=(xano2,-0.08), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )

ax.annotate("$E < E_f$", xy=(xano1,-0.08), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )
plt.show()
