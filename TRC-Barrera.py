##################################################################
#Fisica de Semiconductores 
#Jeison Ivan Roa Mora 
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
figure = Figure() #Figura
canvas = FigureCanvas(figure) #Canvas
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.25)
Vo=1  #Variable para el valor del potencial
a=1.5*10**-9 #Ancho del potencial en metros
En=np.linspace(0.04,3,1000) #Valores de energia de la particula en eV:
T=[]
R=[]
##################################################################
#Declaracion funcion
def den(Vo,a,E,z):
    k=np.sqrt(2*sc.m_e*abs(Vo-E)*sc.e)/sc.hbar
    x=(Vo**2)/(4.0*E)
    if   z==0:
        return x*(np.sinh(k*a)**2)/(Vo-E)
    elif z==1:
        return x*(np.sin(k*a)**2)/(E-Vo)
    else:
        return (sc.m_e*(a**2))*Vo*sc.e/(2*(sc.hbar)**2)
for E in En:
    if E<=Vo:
        x=1/(1+den(Vo,a,E,0))
    elif E>Vo:
        x=1/(1+den(Vo,a,E,1))
    else:
        x=1/(1+den(Vo,a,E,2))
    T.append(x)
    R.append(1-x)
##################################################################
#PLots
plt.title("Coeficientes de Transmision y Reflexion para Barrera de potencial", ha="center", fontsize=15,color="k", fontweight = 'bold')
plt.ylabel(r'E/Vo', fontsize=16,color='k', fontweight = 'bold')
a0 = 1
f0 = 1.5*10**-9
l, = plt.plot(En,T, lw=2, color='crimson')
ll, = plt.plot(En,R, lw=2, color='k')
dot, = plt.plot(En[300],T[300],'bo', markersize=18)
dot2, = plt.plot(En[300],R[300],'ro', markersize=18)
plt.axis([0, 3, -0.3, 1.3])
axcolor = 'gold'
axfreq = plt.axes([0.1, 0.1, 0.65, 0.03], axisbg=axcolor)
axamp = plt.axes([0.1, 0.15, 0.65, 0.03], axisbg=axcolor)
axpunto=plt.axes([0.1, 0.05, 0.65, 0.03], axisbg=axcolor)
aaa = Slider(axfreq, 'Ancho', 1*10**-9, 2*10**-9, valinit=f0)
samp = Slider(axamp, 'Vo', 0.01, 3.0, valinit=a0)
punto = Slider(axpunto, 'E', 0, 1000, valinit=500)
##################################################################
#Declaracion funcion para refrescar y mostrar la grafica
def update(val):
	a = aaa.val
	Vo=samp.val
	p=punto.val
	T=[]
	R=[]
	for E in En:
		if E<=Vo:
			x=1/(1+den(Vo,a,E,0))
		elif E>Vo:
			x=1/(1+den(Vo,a,E,1))
		else:
			x=1/(1+den(Vo,a,E,2))
		T.append(x)
		R.append(1-x)
	l.set_ydata(T)
	ll.set_ydata(R)
	dot.set_data(En[int(p)],T[int(p)])
	dot2.set_data(En[int(p)],R[int(p)])
	fig.canvas.draw_idle()
aaa.on_changed(update)
samp.on_changed(update)
punto.on_changed(update)
plt.plot(En,T,color=plt.cm.gist_rainbow(50))
plt.plot(En,R,color=plt.cm.gist_rainbow(150))
plt.show()
