###################################################################
#Fisica de Semiconductores
Jeison Ivan Roa Mora 
##################################################################
#Importar Librerias
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants   as sc
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.figure import Figure
from matplotlib import colors
from Tkinter import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
##################################################################
#Declaracion de variables, ejes y constantes
figure = Figure()
canvas = FigureCanvas(figure)
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.25)
Vo=1                                            #Valor del potencial en eV
a=1.5*10**-9                                    #Ancho del potencial en metros
En=np.linspace(0.01,1.5,1000)                   #Valores de energia de la particula en eV:
T=[]
R=[]
##################################################################
#Declaracion funcion
for E in En:
    if E<=Vo:
        T.append(0.0)
        R.append(1.0)
    else:
        x=4*np.sqrt(E**2-E*Vo)/((E**(1/2.0)+np.sqrt(E-Vo))**2)
        T.append(x)
        R.append(1-x)
##################################################################
#PLots
plt.title("Coeficientes de Transmision y Reflexion para Escalon de Potencial", ha="center", fontsize=15,color="k", fontweight = 'bold')
plt.ylabel(r'E/Vo', fontsize=16,color='k', fontweight = 'bold')
a0 = 1
f0 = 1.5*10**-9
l, = plt.plot(En,T, lw=2, color='crimson')
ll, = plt.plot(En,R, lw=2, color='k')
dot, = plt.plot(En[300],T[300],'bo', markersize=18)
dot2, = plt.plot(En[300],R[300],'ro', markersize=18)
plt.axis([0, 1.5, -0.3, 1.3])
axcolor = 'gold'
axamp = plt.axes([0.1, 0.15, 0.65, 0.03], axisbg=axcolor)
axpunto=plt.axes([0.1, 0.05, 0.65, 0.03], axisbg=axcolor)
samp = Slider(axamp, 'Vo', 0.01, 1.5, valinit=a0)
punto = Slider(axpunto, 'E', 0, 1000, valinit=100)
##################################################################
#Declaracion funcion para refrescar y mostrar la grafica
def update(val):
	Vo=samp.val
	p=punto.val
	T=[]
	R=[]
	for E in En:
	    if E<=Vo:
		T.append(0.0)
		R.append(1.0)
	    else:
		x=4*np.sqrt(E**2-E*Vo)/((E**(1/2.0)+np.sqrt(E-Vo))**2)
		T.append(x)
		R.append(1-x)
	l.set_ydata(T)
	ll.set_ydata(R)
	dot.set_data(En[int(p)],T[int(p)])
	dot2.set_data(En[int(p)],R[int(p)])
	fig.canvas.draw_idle()
samp.on_changed(update)
punto.on_changed(update)
plt.plot(En,T,color=plt.cm.gist_rainbow(50))
plt.plot(En,R,color=plt.cm.gist_rainbow(150))
plt.show()
