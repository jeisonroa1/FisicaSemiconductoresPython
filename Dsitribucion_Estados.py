## Distribución de estados
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

fig = plt.figure() 
ax = plt.subplot('111', axisbg='w')


#Constantes utilizadas
k  = sc.k
e = sc.e
hb = sc.hbar
m = sc.m_e
a = 6e-15

Cte=(np.sqrt(2.0)*m**(3.0/2)*(a**3))/((np.pi**2)*(hb**3))
print Cte

#Energia en eV relativa
Eo = 5.0; Emin=0; Emax=10.0;
E = np.linspace(Emin,Emax,1000)

ymin=-0.2; ymax=2.7; xval=0.05; yval =2.3;

line, = ax.plot([], [])

plt.ylabel("$\mathcal{N}(E)$")
plt.xlabel("$Energia \ \ (E) \ \ [e\cdot V ]$")
plt.title('$DENSIDAD \ \ DE \ \ ESTADOS \ \ \mathcal{N}(E)$')

def fun(Temp, color1, etq):
    #Funcion de distribucion de Fermi-Dirac
    ax.plot(E,(E**(0.5))*((1)/((np.exp((E-Eo)/(k/e*Temp)))+(1.0))),color=plt.cm.gist_rainbow(color1), label=etq)

    
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

fun(T1,col1, r"$\mathcal{N}(E)_{T = %d K}$" %T1)
fun(T2,col2, r"$\mathcal{N}(E)_{T = %d K}$" %T2)
fun(T3,col3, r"$\mathcal{N}(E)_{T = %d K}$" %T3)
fun(T4,col4, r"$\mathcal{N}(E)_{T = %d K}$" %T4)
fun(T5,col5, r"$\mathcal{N}(E)_{T = %d K}$" %T5)
fun(T6,col6, r"$\mathcal{N}(E)_{T = %d K}$" %T6)

##legend = plt.legend(loc=1, fontsize=12.5)
plt.axis([Emin,Emax,ymin,ymax])

def animate(i): 
    y=(E**(0.5))*((1)/((np.exp((E-Eo)/(k/e*(i*50))))+(1.0)))
    line, = ax.plot(E, y, color=plt.cm.gist_rainbow(i))
    return line,

ani = animation.FuncAnimation(fig, animate,  frames=200, interval=50, blit = True, repeat=True)

legend = plt.legend(loc=1, prop={'size':16})
xano1=((Eo-Emin)/2.0)-0.1
xano2= Eo+((Emax-Eo)/2.0)-0.1
            
ax.annotate(r"$\mathcal{N}(E) = Cte E^{\frac{1}{2}} \frac{1}{e^\frac{E - E_{f}}{kT}+1}$", xy=(0.3,2.25), color='k',
            fontsize=22,
	    horizontalalignment='left',
            verticalalignment='up',
            )
            
ax.annotate(r"$\mathcal{N}(E)_{T = 0K, E<E_{f}} \ \ \propto \ \ E^{\frac{1}{2}}$", xy=(6.5,0.8), color='k',
            fontsize=15,
	    horizontalalignment='left',
            verticalalignment='up',
            )
ax.annotate("$\mathcal{N}(E)_{T = 0K, E>E_{f}} = 0$", xy=(6.5,0.6), color='k',
            fontsize=15,
	    horizontalalignment='left',
            verticalalignment='up',
            )
            
ax.annotate("$E_f$", xy=(Eo-0.2,-0.15), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )

ax.annotate("$E > E_f$", xy=(xano2,-0.15), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )

ax.annotate("$E < E_f$", xy=(xano1,-0.15), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )
plt.show()
