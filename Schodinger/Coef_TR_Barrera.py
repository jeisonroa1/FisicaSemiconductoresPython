## Barrera de Potencial
## Jeison Ivan Roa M.

import numpy 	         as np
import scipy.constants   as sc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.figure import Figure
from matplotlib.widgets import RadioButtons, Slider

from matplotlib import colors
from Tkinter import *
 
#Valor del potencial en eV:
Vo = 1

xmin=0.001; xmax=Vo+2; ymin=-0.5; ymax=2.5; xval=0.05; yval =2.3;

#Ancho de la barrera de potencial en metros (nm):
a=1.5*10**-9

#Valores de energia de la particula en eV:
En0=np.linspace(0.01,Vo-0.001,1000)
En1=np.linspace(Vo+0.001,Vo+2,1000)
x0=[Vo,Vo]
y0=[-0.2,1.2]

def fx1(En0):
    k=np.sqrt(2*sc.m_e*abs(Vo-En0)*sc.e)/sc.hbar
    c=(Vo**2)/(4.0*En0)
    return 1/(1+c*(np.sinh(k*a)**2)/(Vo-En0))

def fx2(En1):
    k=np.sqrt(2*sc.m_e*abs(Vo-En1)*sc.e)/sc.hbar
    c=(Vo**2)/(4.0*En1)
    return 1/(1+c*(np.sin(k*a)**2)/(En1-Vo))

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.16, right=0.9, top=0.75)

rax = plt.axes([0.1, 0.82, 0.8, 0.15], axisbg='white')
radio = RadioButtons(rax, ('Coeficiente de Transmicion (T)', 'Coeficiente de Reflexion (R)',
                           'Coeficiente de Transmicion (T) y Reflexion (R)'))

axdot = plt.axes([0.1, 0.03, 0.8, 0.04], axisbg='white')
sdot = Slider(axdot, "$E \\ [e\cdot V ]$", 0.001, xmax, valinit=0.001)

plt.subplot(111)
plt.xlabel("$Energia \ \ (E) \ \ [e\cdot V ]$")
plt.ylabel(r'$Coeficiente \ \ de \ \ Transmision  \ \ (T)$')
plt.title("$Coeficiente \ \ de \ \ Transmision  \ \ (T) \ \ para \ \ una \ \ barrera \ \ de \ \ potencial  $")
plt.axis([xmin,xmax,ymin,ymax])

#Grafica Coeficiente de transmicion
t0, = plt.plot(En0, fx1(En0), lw=1, color='blue', label= r"$T_{E < V_o} =\frac{4 E (V_o-E)}{4 E (V_o-E) + V_o^2 sinh^2(k_1 a)}$")
t1, = plt.plot(En1, fx2(En1), lw=1, color='blue', label= r"$T_{E > V_o} =\frac{4 E (E-V_o)}{4 E (E-V_o) + V_o^2 sin^2(k_1 a)}$")

legend = plt.legend(loc=1, prop={'size':12.5})

puntox = ax.text(xval,yval, 'Nivel de Energia (E) [eV] = %f' % 0)
puntot = ax.text(xval,yval-0.2, 'Coef. de Transmicion (T) = %f'% 0)
##puntor = ax.text(1.7,-0.1, 'Coef. de Reflexion = %f'% 0.001)

##Punto Coef. de Transmicion
pt, = ax.plot(xmin, fx1(xmin), 'o',color='blue')

#Barrera
b1, = plt.plot(x0, y0, lw=1, c='black', ls='--')

xano1=((Vo-xmin)/2.0)-0.1
xano2= Vo+((xmax-Vo)/2.0)-0.1

ax.annotate("$V_o$", xy=(Vo-0.04,-0.4), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )

ax.annotate("$E > V_o$", xy=(xano2,-0.4), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )

ax.annotate("$E < V_o$", xy=(xano1,-0.4), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )

band_g = 0

def sel_graf(label):
    global band_g, pt, pr
    global puntox, puntot, puntor, band, xmin,xmax,ymin,ymax
    sdot.reset()
    plt.cla()
    plt.axis([xmin,xmax,ymin,ymax])

    #Barrera
    b1, = plt.plot(x0, y0, lw=1, c='black', ls='--')

    puntox = ax.text(xval,yval, 'Nivel de Energia (E) [eV] = %f' %0)

    sel_coef = {'Coeficiente de Transmicion (T)': 0,
              'Coeficiente de Reflexion (R)': 1,
              'Coeficiente de Transmicion (T) y Reflexion (R)': 2}
    
    label1 = sel_coef[label]
    
    if label1 == 0:
        band_g = 0

        #Grafica Coeficiente de transmicion
        t0, = plt.plot(En0, fx1(En0), lw=1, color='blue', label= r"$T_{E < V_o} =\frac{4 E (V_o-E)}{4 E (V_o-E) + V_o^2 sinh^2(k_1 a)}$")
        t1, = plt.plot(En1, fx2(En1), lw=1, color='blue', label= r"$T_{E > V_o} =\frac{4 E (E-V_o)}{4 E (E-V_o) + V_o^2 sin^2(k_1 a)}$")

        #Grafica Punto Coef. de Transmicion
        pt, = ax.plot(xmin, fx1(xmin), 'o',color='blue')
        
        plt.xlabel("$Energia \ \ (E) \ \ [e\cdot V ]$")
        plt.ylabel("$Coeficiente \ \ de \ \ Transmision  \ \ (T)$")
        plt.title("$Coeficiente \ \ de \ \ Transmision  \ \ (T) \ \ para \ \ una \ \ barrera \ \ de \ \ potencial  $")

        puntot = ax.text(xval,yval-0.2, 'Coef. de Transmicion (T) = %f'% 0)
        
    elif label1 == 1:
        band_g = 1

        #Grafica Coeficiente de Reflexion
        r0, = plt.plot(En0, 1-fx1(En0), lw=1, color='red', label=r"$R_{E < V_o} =1-\frac{4 E (V_o-E)}{4 E (V_o-E) + V_o^2 sinh^2(k_1 a)}$")
        r1, = plt.plot(En1, 1-fx2(En1), lw=1, color='red', label=r"$R_{E > V_o} =1-\frac{4 E (E-V_o)}{4 E (E-V_o) + V_o^2 sin^2(k_1 a)}$")

        #Grafica Punto Coef. de Reflexion
        pr, = ax.plot(xmin, 1-fx1(xmin), 'o',color='red')
        
        plt.xlabel("$Energia \ \ (E) \ \ [e\cdot V ]$")
        plt.ylabel("$Coeficiente \ \ de \ \ Reflexion  \ \ (R)$")
        plt.title("$Coeficiente \ \ de \ \ Reflexion \ \ (R) \ \ para \ \ una \ \ barrera \ \ de \ \ potencial  $")

        puntor = ax.text(xval,yval-0.2, 'Coef. de Reflexion (R)= %f'% 1)

    elif label1 == 2:
        band_g = 2

        #Grafica Punto Coef. de Transmicion
        pt, = ax.plot(xmin, fx1(xmin), 'o',color='blue')
        #Grafica Punto Coef. de Reflexion
        pr, = ax.plot(xmin, 1-fx1(xmin), 'o',color='red')

        #Grafica Coeficiente de Reflexion
        t0, = plt.plot(En0, fx1(En0), lw=1, color='blue', label= r"$T_{E < V_o} =\frac{4 E (V_o-E)}{4 E (V_o-E) + V_o^2 sinh^2(k_1 a)}$")
        t1, = plt.plot(En1, fx2(En1), lw=1, color='blue', label= r"$T_{E > V_o} =\frac{4 E (E-V_o)}{4 E (E-V_o) + V_o^2 sin^2(k_1 a)}$")
        r0, = plt.plot(En0, 1-fx1(En0), lw=1, color='red', label= r"$R_{E < V_o} =1-T_{E < V_o}$")
        r1, = plt.plot(En1, 1-fx2(En1), lw=1, color='red', label= r"$R_{E > V_o} =1-T_{E > V_o}$")

        puntot = ax.text(xval,yval-0.2, 'Coef. de Transmicion (T) = %f'% 0)
        puntor = ax.text(xval,yval-0.4, 'Coef. de Reflexion (R)= %f'% 1)

        plt.xlabel("$Energia \ \ (E) \ \ [e\cdot V ]$")
        plt.ylabel("$Coeficiente \ \ de \ \ Transmision \ \ (T) \ \ y \ \ Reflexion \ \ (R)$")
        plt.title("$Coeficiente \ \ de \ \ Transmision \ \ (T) \ \ y \ \ Reflexion \ \ (R) \ \ para \ \
                    una \ \ barrera \ \ de \ \ potencial $")

    ax.annotate("$V_o$", xy=(Vo-0.04,-0.4), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )

    ax.annotate("$E > V_o$", xy=(xano2,-0.4), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )

    ax.annotate("$E < V_o$", xy=(xano1,-0.4), color='black',
	    horizontalalignment='left',
            verticalalignment='up',
            )
    
    legend = plt.legend(loc=1, prop={'size':12.5})
    
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
    plt.draw()
    print band_g
    
radio.on_clicked(sel_graf)

def graf_p(val):
    global puntox, puntot, puntor, band
    dot = sdot.val

    if band_g==0:  
        if dot<Vo:
            pt.set_xdata(dot)
            pt.set_ydata(fx1(dot))
            t = fx1(dot)            
        elif dot>Vo:
            pt.set_xdata(dot)
            pt.set_ydata(fx2(dot))
            t = fx2(dot)
        puntox.remove()
        puntot.remove()
        puntox = ax.text(xval,yval, 'Nivel de Energia (E) [eV] = %f' % dot)
        puntot = ax.text(xval,yval-0.2, 'Coef. de Transmicion (T) = %f'% t)
        
    elif band_g==1:
        if dot<Vo:
            pr.set_xdata(dot)
            pr.set_ydata(1-fx1(dot))
            r = 1-fx1(dot)
        elif dot>Vo:
            pr.set_xdata(dot)
            pr.set_ydata(1-fx2(dot))
            r = 1-fx2(dot)
        puntox.remove()
        puntor.remove()
        puntox = ax.text(xval,yval, 'Nivel de Energia (E) [eV] = %f' % dot)
        puntor = ax.text(xval,yval-0.2, 'Coef. de Reflexion (R) = %f'% r)
            
    elif band_g==2:
        if dot<Vo:
            pt.set_xdata(dot)
            pt.set_ydata(fx1(dot))
            pr.set_xdata(dot)
            pr.set_ydata(1-fx1(dot))
            t = fx1(dot)
            r = 1-fx1(dot)
        elif dot>Vo:
            pt.set_xdata(dot)
            pt.set_ydata(fx2(dot))
            pr.set_xdata(dot)
            pr.set_ydata(1-fx2(dot))
            t = fx2(dot)
            r = 1-fx2(dot)
        puntox.remove()
        puntot.remove()
        puntor.remove()
        puntox = ax.text(xval,yval, 'Nivel de Energia (E) [eV] = %f' % dot)
        puntot = ax.text(xval,yval-0.2, 'Coef. de Transmicion (T) = %f'% t)
        puntor = ax.text(xval,yval-0.4, 'Coef. de Reflexion (R)= %f'% r)
            
    else:
        pass
    
    fig.canvas.draw()   
sdot.on_changed(graf_p)

plt.show()
