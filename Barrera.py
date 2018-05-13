###################################################################
#Fisica de Semiconductores - Universidad Distrital Francisco Jose de Caldas
#Jeison Ivan Roa Mora - 20131005056
#Gustavo Rivas Gutierrez - 20132005036 

###################################################################
#Se importan librerias
import matplotlib.animation as animation
import numpy as np    
from scipy import signal as sp
import matplotlib.pyplot as plt         
###################################################################
#Definicion de variables
Vo=16 #Altura del escalon
m=1
E=10
D=10
a=2.4
A=5
B=10
C=1.25
D=2.5
x = np.arange(-1, 5, 0.001)
K1=np.sqrt((2*m*E))
K2=np.sqrt((2*m*(Vo-E)))
###################################################################        
#Definicion del rango de trabajo, e inicializacion de datos
def init():
    ax.set_ylim(-25, 25)
    ax.set_xlim(-4,10)
    del xdata[:]
    del ydata[:]
    line.set_data(xdata, ydata)
    return line,
###################################################################
#Funcion refrescar datos
def data_gen(t=-4):
    cnt = 0
    while cnt < 1000:
        cnt += 1
        t += 0.07
        if t<=0:
            yield t, A*np.exp(1j*t*K1)+B*np.exp(-1j*K1*t)
        elif t<a:
            yield t, D*np.exp(-K2*t)
        else:
            yield t, C*np.exp(1j*K1*t)+D*np.exp(-1j*K1*t)
###################################################################
#Calibracion de eje
fig, ax = plt.subplots()
line, = plt.plot([], [], lw=2, color='blue')
ax.grid()
xdata, ydata = [], []
###################################################################
#Graficar y actualizar datos
def run(data):
    plt.title("Grafica Potencial de barrera",fontsize=16,fontweight='bold', color='k')
    plt.xlabel("$x$",fontsize=20)
    plt.ylabel("$|\psi(x)|$",fontsize=20)
    t, y = data
    xdata.append(t)
    ydata.append(y)
    xmin, xmax = ax.get_xlim()
    plt.plot([-4,0],[0,0])
    plt.plot([0,0],[0,Vo])
    plt.plot([0,2.4],[Vo,Vo])
    plt.plot([2.4,2.4],[0,Vo])
    plt.plot([2.4,10],[0,0])
    line.set_data(xdata, ydata)    
    return line,
###################################################################
#Refrescar pantalla animada
ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10,repeat=False, init_func=init)
plt.show()

