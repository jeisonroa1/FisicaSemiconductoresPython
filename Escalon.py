###################################################################
#Fisica de Semiconductores - Universidad Distrital Francisco Jose de Caldas
#Jeison Ivan Roa Mora - 20131005056
#Gustavo Rivas Gutierrez - 20132005036 

###################################################################
# Importar librerias
import numpy as np                       
import time 
import matplotlib.pyplot as plt          
import scipy.constants as sc                
from matplotlib.widgets import Slider       
from matplotlib.widgets import Button      
import matplotlib.animation as animation    
import matplotlib.gridspec as gridspec      
###################################################################
# Defincion de constantes
global Psi
V0 = 5*sc.eV                                # Voltaje de Escalon
E = np.linspace(-2*V0,2*V0,1e3)             # Energia de la particula
x = np.linspace(-10e-10,10e-10,1e3)         # Posicion de la particula
t = np.linspace(0,2*np.pi,50)               # Tiempo
x1 = np.linspace(-10e-10,0,500)             # X1
x2 = np.linspace(0,10e-10,500)              # X2

###################################################################
# Definicion funcion para X<0
def Pot(x):                                 
	if (x < 0):
		return 0
	else:
		return V0
###################################################################
# Definicion funcion de onda 
def var_amp1(A): 
	global Psi                                     

	if (V0>(A*sc.eV)):
		k2 = np.sqrt(2*sc.m_e*(V0 - A*sc.eV))/sc.hbar
		k3 = k2
	else:
		k3 = np.sqrt(2*sc.m_e*(A*sc.eV - V0))/sc.hbar
		k2=k3
	k1 = np.sqrt(2*sc.m_e*A*sc.eV)/sc.hbar
	if(V0 > A*sc.eV):
		psi1 = np.cos(k1*x1) - (k2/k1)*np.sin(k1*x1)
		psi2 = np.exp(-k2*x2)
	else:
		psi1 = np.exp(1j*k1*x1) + (k1-k3)*(np.exp(-1j*k1*x1))/(k1+k3)
		psi2 = (2*k1)*np.exp(1j*k3*x2)/(k1+k3)
	
	Psi = np.concatenate((psi1,psi2))
	imsm = [];
	for add in t:
		z = Psi * np.exp(-1j*add)
		mag = z * np.conj(z);
		imsm.append(plt.plot(x,mag,color = 'tomato'))
	return imsm
###################################################################
#Parte real de la funcion de onda
def var_amp2(A):
	global Psi    
	imsr = [];
	for add in t:
		z = Psi * np.exp(-1j*add)
		imsr.append(plt.plot(x,np.real(z),color = 'lime'))
	return imsr
###################################################################
#Parte Imaginaria de la funcion de onda
def var_amp3(A):
	global Psi        
	imsi = [];
	for add in t:
		z = Psi * np.exp(-1j*add)
		imsi.append(plt.plot(x,np.imag(z),color = 'gold'))
	return imsi
###################################################################
#Funcion para el Slider de la barra de energia
def update(val):                        
	energy = barra.val
	aux.set_xdata(energy)
	f4.set_text('Energia = %e eV'%energy)
 
	if energy*sc.eV > V0:
		f3.set_text(r'$V_0$ = %e eV'%Pot(1))
	else:
		f3.set_text(r'$V_0$ = %e eV'%Pot(-1))
	
	global im_ani1,im_ani2,im_ani3,imsm,imsr,imsi     
	grafica_sim1()
	if (energy != 0):
		imsm = var_amp1(energy)
		im_ani1 = animation.ArtistAnimation(anim,imsm,interval=1,repeat_delay=1,blit=False)
	 
	if (energy != 0):
		imsr = var_amp2(energy)
		im_ani2 = animation.ArtistAnimation(anim,imsr,interval=1,repeat_delay=1,blit=False)      
		
	if (energy != 0):
		imsi = var_amp3(energy)
		im_ani3 = animation.ArtistAnimation(anim,imsi,interval=1,repeat_delay=1,blit=False)  

	plt.draw()
	time.sleep(0.1)
###################################################################
#Animacion de las graficas
def grafica_sim1():                          
	plt.subplot(G[0,:])
	plt.cla()    
	plt.plot(x,Potencial,'y',linewidth=2,color='k')
	plt.text(5e-10,-1.4,'______________',color='k',fontsize=16)
	plt.text(5e-10,-3.1,'______________',color='k',fontsize=16)
	plt.text(5e-10,-1.9,'|',color='k',fontsize=16)
	plt.text(5e-10,-2.3,'|',color='k',fontsize=16)
	plt.text(5e-10,-2.7,'|',color='k',fontsize=16)
	plt.text(5e-10,-3.1,'|',color='k',fontsize=16)
	plt.text(7.35e-10,-1.9,'|',color='k',fontsize=16)
	plt.text(7.35e-10,-2.3,'|',color='k',fontsize=16)
	plt.text(7.35e-10,-2.7,'|',color='k',fontsize=16)
	plt.text(7.35e-10,-3.1,'|',color='k',fontsize=16)
	plt.text(5e-10,-2,r' Mag($\Psi$(x,t))',color='tomato',fontsize=16)
	plt.text(5e-10,-2.5,r' Real($\Psi$(x,t))',color='lime',fontsize=16)
	plt.text(5e-10,-3,r' Imag($\Psi$(x,t))',color='gold',fontsize=16)
	plt.title("Grafica Dinamica - Escalon de Potencial", ha="center", fontsize=24,color="k", fontweight = 'bold')
	plt.ylabel(r'$\Psi$(x,t)', fontsize=20,color='k', fontweight = 'bold')
	plt.ylim([-5,8]); 
	plt.grid(True)
###################################################################
# Funciones a graficar
Potencial = []
for i in x:
	Potencial.append(Pot(i)/sc.eV)
anim = plt.figure(figsize=[20,8],facecolor='w',frameon=False)
G = gridspec.GridSpec(1,1)
###################################################################
# Grafica Estatica
aux ,= plt.plot(np.zeros(np.size(E)),np.linspace(-6,6,np.size(E)),'--',color='black',linewidth=2)
s = 'Enegia (eV)'
plt.ylim([-0.5,1.5]); plt.xlim([0,10])
plt.grid(True)
f1 = plt.text(4,-0.2,'',fontsize = 12, color = 'r')
f2 = plt.text(4,-0.3,'',fontsize = 12, color = 'b')
f3 = plt.text(4,-0.4,'',fontsize = 12, color = 'y')
f4 = plt.text(4,-0.1,'',fontsize = 12, color = 'k')
###################################################################
# Graficas animadas
grafica_sim1()
###################################################################
# Definicion y metodo barra interactiva y boton
barra = Slider(plt.axes([0.2,0.015,0.55,0.03], axisbg='darkolivegreen'), r'Energia (eV)',0,5,valinit=0,color='lime')
barra.on_changed(update)

plt.show()
