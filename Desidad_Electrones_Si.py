## Densidad electrones Si
## Jeison Ivan Roa M.


import numpy             as np
import matplotlib.pyplot as plt
import scipy.constants   as sc
from matplotlib.widgets import RadioButtons, Slider


#--------------------------------------------------------------------------------------------------------
#DEL APENDICE DEL SZE:

semic=['Ge','Si']
cl=['b','k','mediumblue']
Ego=[0.7437,1.170,1.519,1.421,2.34,3.47]
a=[4.774,4.73,5.405,4.9,6,7.7]
b=[235,636,204,327,460,600]

#Gap a 300K:
Eg=[0.66,1.12,1.42,1.35,2.26,3.36]

#Densidad de portadores intrinsecos a 300K:
ni=[2.4*10**13,1.45*10**10,1.79*10**6, 1.3*10**7,2*10**6,1.0*10**10]

#Concentracion de donadores (cm^(-3))
ND=1e16

#Concentracion de aceptadores (cm^(-3))
NA=1e10

gD=2
ED=0.05
EF=1.07
mn=1.08
m0=9.1e-31

#--------------------------------------------------------------------------------------------------------
#Densidad de electrones n (cm^(-3))
def n_Si(T):
    sqrtNcNv_Si = ni[1]*np.exp(Eg[1]*sc.e/(2*sc.k*300))
    mn_mp_Si = ((sqrtNcNv_Si)/2.0)*(((sc.h**2)/(2*sc.pi*sc.k*300))**(3/2.0))
    Egg_Si=Ego[1]-((a[1]*10**-4)*T**2)/(T+b[1])
    ni_Si=2*((2*sc.pi*sc.k/(sc.h**2))**(3/2.0))*(T**(3/2.0))*mn_mp_Si*np.exp(-Egg_Si*sc.e/(2*T*sc.k))
    return 0.5*(ND-NA+np.sqrt((ND-NA)**2+4*ni_Si**2))

#--------------------------------------------------------------------------------------------------------
#Densidad de electrones n en la region extrinseca y de congelamiento (Freeze-out and extrinsic region) (cm^(-3))
def n_for(T_for):
    Nc=(2*2*np.pi*mn*m0*sc.k*T_for/sc.h**2)**(3.0/2)*0.000001
##    Nc3=(2*2*np.pi*mn*m0*sc.k*300/sc.h**2)**(3.0/2)*0.000001
##    print Nc3
    Nast=Nc/2*np.exp(-ED*sc.e/sc.k/T_for)
    return -(Nast+NA)/2+np.sqrt(((Nast+NA)/2)**2+(ND-NA)*Nast)

#Temperatura de caambio de region T_cr
sqrtNcNv_Si = ni[1]*np.exp(Eg[1]*sc.e/(2*sc.k*300))
T_cr=Eg[1]*sc.e/(sc.k*np.log(sqrtNcNv_Si**2)-2*sc.k*np.log(ND-NA))
print T_cr

#--------------------------------------------------------------------------------------------------------
#Usando la variacion de Eg(T):

#Densidad de portadores intrinsecos en funcion de T:  
def Ni_Si(T):
    #Calculo de masas efectivas del semiconductor a partir de ni(300k)
    sqrtNcNv_Si = ni[1]*np.exp(Eg[1]*sc.e/(2*sc.k*300))
    mn_mp_Si = ((sqrtNcNv_Si)/2.0)*(((sc.h**2)/(2*sc.pi*sc.k*300))**(3/2.0))
    Egg_Si=Ego[1]-((a[1]*10**-4)*T**2)/(T+b[1])
    return 2*((2*sc.pi*sc.k/(sc.h**2))**(3/2.0))*(T**(3/2.0))*mn_mp_Si*np.exp(-Egg_Si*sc.e/(2*T*sc.k))

#-------------------------------------------------------------------------------------------------------
#Rangos de Temperaturas
Tin = 300
T0 = 0.000001
Tf = 700
Tint = T_cr =400
Tint = 400
ft = 15
T = np.linspace(Tint,Tf,1000)
T_for = np.linspace(T0,Tint,1000)

#-------------------------------------------------------------------------------------------------------
#Graficas
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9)

xmin=T0; xmax=Tf; ymin=0; ymax=3e16; xval=10; yval =3e16;
axdot = plt.axes([0.1, 0.03, 0.8, 0.04], axisbg='w')
sdot = Slider(axdot, "$T \\ [K]$", xmin, xmax, valinit=Tin)

plt.subplot(111)

ax.plot(T,n_Si(T), lw=1, color = cl[0], label=r"$n=\frac{1}{2}\left[N_D-N_A+\sqrt{(N_D-N_A)^2+4*n_i^2}\right]$")
ax.plot(T,Ni_Si(T), ls='--', lw=1, color = cl[1], label=r"$n \ \ = \ \ n_i \ \ , \ \ N_D=N_A=10^6 cm^{-3}$")
ax.plot(T_for,n_for(T_for), ls='--', lw=1, color = cl[2],
        label=r"$n_{\ \ Freeze-out \ \ and \ \ extrinsic \ \ region}=-\frac{N^*+N_A}{2}+\sqrt{\frac{(N^*+N_A)^2}{4}+N^*(N_D-N_A)}, \ \ donde \ \ N^*=\frac{N_C}{2}e^{-\frac{E_C-E_D}{kT}}$")


ax1 = ax.twiny()
ax.set_xlabel(r"$T \ \ (K)$")
ax1.set_xlabel(r"$T \ \ (^{o}C)$")
ax1.set_xlim(0, 600)
ax1.set_xticks([73, 173, 273, 373, 473, 573, 673])
ax1.set_xticklabels(['-200','-100','0','100','200','300','400'])

#Anotaciones
F_reg = ax.text(15,yval-1.9e16, r"$Freeze-out \ \ region$", fontsize=ft)
Ext_reg = ax.text(350,yval-1.9e16, r"$Extrinsic \ \ region$", fontsize=ft)
Int_reg = ax.text(570,yval-1e16, r"$Intrinsic \ \ region$", fontsize=ft)
    
##import matplotlib.transforms as mtransforms
##trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
##ax.fill_between(T, 0, 1, where = T > T_cr, facecolor='blue', alpha=0.05, transform=trans)
##ax.fill_between(T_for, 0, 1, where = T_for < 140, facecolor='red', alpha=0.05, transform=trans)

#Puntos
p1, = ax.plot(300, n_for(300), 'o', color = cl[0])

def graf_p(val):
    global puntox, punto1, punto2, punto3, punto4, punto5, punto6
    dot = sdot.val
    dot = round(sdot.val,2)
    if dot<T_cr:
        p1.set_xdata(dot)
        p1.set_ydata(n_for(dot))
        puntox.remove()
        punto1.remove()
        punto1 = ax.text(xval, yval-1.2e16, r"$Densidad \ \ de \ \ Electrones \ \ n \ \ [cm^{-3}] \ \ = \ \ %.3e$" % n_for(dot), color= 'b', fontsize=ft) 
##    print dot
    else:
        p1.set_xdata(dot)
        p1.set_ydata(n_Si(dot))
        puntox.remove()
        punto1.remove()
        punto1 = ax.text(xval, yval-1.2e16, r"$Densidad \ \ de \ \ Electrones \ \ n \ \ [cm^{-3}] \ \ = \ \ %.3e$" % n_Si(dot), color= 'b', fontsize=ft) 
    puntox = ax.text(xval, yval-1e16, r"$Temperatura \ \ (T) \ \ [K] \ \ = \ \ %.2f$" % dot, fontsize=ft)
    fig.canvas.draw()   
sdot.on_changed(graf_p)

puntox = ax.text(xval, yval-1e16, r"$Temperatura \ \ (T) \ \ [K] \ \ = \ \ %.2f$" % Tin, fontsize=ft)
punto1 = ax.text(xval, yval-1.2e16, r"$Densidad \ \ de \ \ Electrones \ \ n \ \ [cm^{-3}] \ \ = \ \ %.3e$" % n_for(Tin), color= 'b', fontsize=ft)

ax.legend(loc=2, prop={'size':15})
ax.set_ylabel(r"$Densidad \ \ de \ \ Electrones \ \ n \ \ (cm^{-3}), \ \ Si$",{'fontsize': 20})
plt.axis([xmin, xmax, ymin, ymax])
plt.show()
