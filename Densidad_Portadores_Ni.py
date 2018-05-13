## Densidad de portadores
## Jeison Ivan Roa M.

import numpy             as np
import matplotlib.pyplot as plt
import scipy.constants   as sc
from matplotlib.widgets import RadioButtons, Slider

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9)

xmin=0; xmax=2000; ymin=1; ymax=1e22; xval=50; yval =3.1;

Tamb=[300,300]
yTamb=[ymin, ymax]

Tin = 300
ft=12

axdot = plt.axes([0.1, 0.03, 0.8, 0.04], axisbg='w')
sdot = Slider(axdot, "$T \\ [K]$", 0, xmax, valinit=Tin)

plt.subplot(111)

semic=['Ge','Si','GaAs','InP','GaP','GaN']
color=['b','r','g','m','y','k']
cl=['b','r','g','m','y','k']
Ego=[0.7437,1.170,1.519,1.421,2.34,3.47]
a=[4.774,4.73,5.405,4.9,6,7.7]
b=[235,636,204,327,460,600]

#Del apendice del sze:

#Gap a 300K:
Eg=[0.66,1.12,1.42,1.35,2.26,3.36]

#Densidad de portadores intrinsecos a 300K:
ni=[2.4*10**13,1.45*10**10,1.79*10**6, 1.3*10**7,2*10**6,1.0*10**10]
T=np.linspace(100,2000)

#Calculo de masas efectivas del semiconductor a partir de ni(300k)
#Usando la variacion de Eg(1000/T):

#Densidad de portadores intrinsecos en funcion de 1000/T:  
def Ni_Ge(T):
    sqrtNcNv_Ge = ni[0]*np.exp(Eg[0]*sc.e/(2*sc.k*300))
    mn_mp_Ge = ((sqrtNcNv_Ge)/2.0)*(((sc.h**2)/(2*sc.pi*sc.k*300))**(3/2.0))
    Egg_Ge=Ego[0]-((a[0]*10**-4)*T**2)/(T+b[0])
    return 2*((2*sc.pi*sc.k/(sc.h**2))**(3/2.0))*(T**(3/2.0))*mn_mp_Ge*np.exp(-Egg_Ge*sc.e/(2*T*sc.k))

def Ni_Si(T):
    sqrtNcNv_Si = ni[1]*np.exp(Eg[1]*sc.e/(2*sc.k*300))
    mn_mp_Si = ((sqrtNcNv_Si)/2.0)*(((sc.h**2)/(2*sc.pi*sc.k*300))**(3/2.0))
    Egg_Si=Ego[1]-((a[1]*10**-4)*T**2)/(T+b[1])
    return 2*((2*sc.pi*sc.k/(sc.h**2))**(3/2.0))*(T**(3/2.0))*mn_mp_Si*np.exp(-Egg_Si*sc.e/(2*T*sc.k))

def Ni_GaAs(T):
    sqrtNcNv_GaAs = ni[2]*np.exp(Eg[2]*sc.e/(2*sc.k*300))
    mn_mp_GaAs = ((sqrtNcNv_GaAs)/2.0)*(((sc.h**2)/(2*sc.pi*sc.k*300))**(3/2.0))
    Egg_GaAs=Ego[2]-((a[2]*10**-4)*T**2)/(T+b[2])
    return 2*((2*sc.pi*sc.k/(sc.h**2))**(3/2.0))*(T**(3/2.0))*mn_mp_GaAs*np.exp(-Egg_GaAs*sc.e/(2*T*sc.k))

def Ni_InP(T):
    sqrtNcNv_InP = ni[2]*np.exp(Eg[2]*sc.e/(2*sc.k*300))
    mn_mp_InP = ((sqrtNcNv_InP)/2.0)*(((sc.h**2)/(2*sc.pi*sc.k*300))**(3/2.0))
    Egg_InP=Ego[3]-((a[3]*10**-4)*T**2)/(T+b[3])
    return 2*((2*sc.pi*sc.k/(sc.h**2))**(3/2.0))*(T**(3/2.0))*mn_mp_InP*np.exp(-Egg_InP*sc.e/(2*T*sc.k))

def Ni_GaP(T):
    sqrtNcNv_GaP = ni[2]*np.exp(Eg[2]*sc.e/(2*sc.k*300))
    mn_mp_GaP = ((sqrtNcNv_GaP)/2.0)*(((sc.h**2)/(2*sc.pi*sc.k*300))**(3/2.0))
    Egg_GaP=Ego[4]-((a[4]*10**-4)*T**2)/(T+b[4])
    return 2*((2*sc.pi*sc.k/(sc.h**2))**(3/2.0))*(T**(3/2.0))*mn_mp_GaP*np.exp(-Egg_GaP*sc.e/(2*T*sc.k))

def Ni_GaN(T):
    sqrtNcNv_GaN = ni[2]*np.exp(Eg[2]*sc.e/(2*sc.k*300))
    mn_mp_GaN = ((sqrtNcNv_GaN)/2.0)*(((sc.h**2)/(2*sc.pi*sc.k*300))**(3/2.0))
    Egg_GaN=Ego[5]-((a[5]*10**-4)*T**2)/(T+b[5])
    return 2*((2*sc.pi*sc.k/(sc.h**2))**(3/2.0))*(T**(3/2.0))*mn_mp_GaN*np.exp(-Egg_GaN*sc.e/(2*T*sc.k))


#Graficas
g1, = plt.semilogy(T,Ni_Ge(T), lw=1, color = cl[0], label=r"$n_{i \ \ Ge, \ \ E_g(T)_{%s}= %.3f-\frac{%.3fT^2}{T+%.0f}}$" %(semic[0],Ego[0],a[0],b[0]))
g2, = plt.semilogy(T,Ni_Si(T), lw=1, color = cl[1], label=r"$n_{i \ \ Si, \ \ E_g(T)_{%s}= %.3f-\frac{%.3fT^2}{T+%.0f}}$" %(semic[1],Ego[1],a[1],b[1]))
g3, = plt.semilogy(T,Ni_GaAs(T), lw=1, color = cl[2], label=r"$n_{i \ \ GaAs, \ \ E_g(T)_{%s}= %.3f-\frac{%.3fT^2}{T+%.0f}}$" %(semic[2],Ego[2],a[2],b[2]))
g4, = plt.semilogy(T,Ni_InP(T), lw=1, color = cl[3], label=r"$n_{i \ \ InP, \ \ E_g(T)_{%s}= %.3f-\frac{%.3fT^2}{T+%.0f}}$" %(semic[3],Ego[3],a[3],b[3]))
g5, = plt.semilogy(T,Ni_GaP(T), lw=1, color = cl[4], label=r"$n_{i \ \ GaP, \ \ E_g(T)_{%s}= %.3f-\frac{%.3fT^2}{T+%.0f}}$" %(semic[4],Ego[4],a[4],b[4]))
g6, = plt.semilogy(T,Ni_GaN(T), lw=1, color = cl[5], label=r"$n_{i \ \ GaN, \ \ E_g(T)_{%s}= %.3f-\frac{%.3fT^2}{T+%.0f}}$" %(semic[5],Ego[5],a[5],b[5]))

#Puntos
p1, = ax.plot(Tin, Ni_Ge(Tin), 'o', color = cl[0])
p2, = ax.plot(300, Ni_Si(300), 'o', color = cl[1])
p3, = ax.plot(300, Ni_GaAs(300), 'o', color = cl[2])
p4, = ax.plot(300, Ni_InP(300), 'o', color = cl[3])
p5, = ax.plot(300, Ni_GaP(300), 'o', color = cl[4])
p6, = ax.plot(300, Ni_GaN(300), 'o', color = cl[5])

#Barrera
##b1, = plt.plot(Tamb, yTamb, lw=1, c='black', ls='--')

def graf_p(val):
    global puntox, punto1, punto2, punto3, punto4, punto5, punto6
    dot = sdot.val
    dot = round(sdot.val,2)
##    print dot 
    p1.set_xdata(dot)
    p1.set_ydata(Ni_Ge(dot))
    p2.set_xdata(dot)
    p2.set_ydata(Ni_Si(dot))
    p3.set_xdata(dot)
    p3.set_ydata(Ni_GaAs(dot))
    p4.set_xdata(dot)
    p4.set_ydata(Ni_InP(dot))
    p5.set_xdata(dot)
    p5.set_ydata(Ni_GaP(dot))
    p6.set_xdata(dot)
    p6.set_ydata(Ni_GaN(dot))
    puntox.remove()
    punto1.remove()
    punto2.remove()
    punto3.remove()
    punto4.remove()
    punto5.remove()
    punto6.remove()
    puntox = ax.text(xval, 1e21, r"$Temperatura \ \ (T) \ \ [K] \ \ = \ \ %.2f$" % dot, fontsize=ft)
    punto1 = ax.text(xval, 1e20, r"$n_{i \ \ Ge} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_Ge(dot), color= 'b', fontsize=ft)
    punto2 = ax.text(xval, 1e19, r"$n_{i \ \ Si} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_Si(dot), color= 'r', fontsize=ft)
    punto3 = ax.text(xval, 1e18, r"$n_{i \ \ GaAs} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_GaAs(dot), color= 'g', fontsize=ft)
    punto4 = ax.text(xval, 1e17, r"$n_{i \ \ InP} \ \ [cm^{-3}] \ \ = \ \ %f.3$" % Ni_InP(dot), color= 'm', fontsize=ft)
    punto5 = ax.text(xval, 1e16, r"$n_{i \ \ GaP} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_GaP(dot), color= 'y', fontsize=ft)
    punto6 = ax.text(xval, 1e15, r"$n_{i \ \ GaN} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_GaN(dot), color= 'k', fontsize=ft)
    fig.canvas.draw()   
sdot.on_changed(graf_p)

puntox = ax.text(xval, 1e21, r"$Temperatura \ \ (T) \ \ [K] \ \ = \ \ %.2f$" % Tin, fontsize=ft)
punto1 = ax.text(xval, 1e20, r"$n_{i \ \ Ge} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_Ge(Tin), color= 'b', fontsize=ft)
punto2 = ax.text(xval, 1e19, r"$n_{i \ \ Si} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_Si(Tin), color= 'r', fontsize=ft)
punto3 = ax.text(xval, 1e18, r"$n_{i \ \ GaAs} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_GaAs(Tin), color= 'g', fontsize=ft)
punto4 = ax.text(xval, 1e17, r"$n_{i \ \ InP} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_InP(Tin), color= 'm', fontsize=ft)
punto5 = ax.text(xval, 1e16, r"$n_{i \ \ GaP} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_GaP(Tin), color= 'y', fontsize=ft)
punto6 = ax.text(xval, 1e15, r"$n_{i \ \ GaN} \ \ [cm^{-3}] \ \ = \ \ %.3f$" % Ni_GaN(Tin), color= 'k', fontsize=ft)

plt.legend = plt.legend(loc=4, prop={'size':20})
plt.title(r"$Densidad \ \ intrinseca \ \ de \ \ portadores \ \ n_i$")
plt.ylabel(r"$Densidad \ \ intrinseca \ \ de \ \ portadores \ \ n_i \ \ (cm^{-3})$")
plt.xlabel(r"$T \ \ (K)$")
##plt.semilogy([300,300],[10**6,10**25],'--')
plt.axis([xmin, xmax, ymin, ymax])
##plt.ylim(1e1,1e25)
##plt.grid()
plt.show()
