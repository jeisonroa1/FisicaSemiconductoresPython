############################################################
#Fisica de semiconductores
#jeison Ivan Mora Roa 20131005056
#Gustavo Rivas Gutierrez 20132005036
###############################################################
import numpy as np
import pylab
pylab.ion()

def Gaussian(x,t,sigma):
    """  A Gaussian curve.
        x = Variable
        t = time shift
        sigma = standard deviation      """
    return np.exp(-(x-t)**2/(2*sigma**2))
def free(npts):
    "Free particle."
    return np.zeros(npts)
def step(npts,v0):
    "Potential step"
    v = free(npts)
    v[npts/2:] = v0
    return v
def barrier(npts,v0,thickness):
    "Barrier potential"
    v = free(npts)
    v[npts/2:npts/2+thickness] = v0
    return v
def fillax(x,y,*args,**kw):
    """Fill the space between an array of y values and the x axis.
    All args/kwargs are passed to the pylab.fill function.
    Returns the value of the pylab.fill() call.
    """
    xx = np.concatenate((x,np.array([x[-1],x[0]],x.dtype)))
    yy = np.concatenate((y,np.zeros(2,y.dtype)))
    return pylab.fill(xx, yy, *args,**kw)
#
N    = 1200     
T    = 5*N     
                
Tp   = 50       
dx   = 1.0e0    
m    = 1.0e0    
hbar = 1.0e0    
X    = dx*np.linspace(0,N,N)       
V0   = 0   
THCK = 0   
POTENTIAL = 'free'

sigma = 40.0 
x0 = round(N/2) - 10*sigma 
k0 = np.pi/20 
E = (hbar**2/2.0/m)*(k0**2+0.5/sigma**2)
################################################################

if POTENTIAL=='free':
    V = free(N)
elif POTENTIAL=='step':
    V = step(N,V0)
elif POTENTIAL=='barrier':
    V = barrier(N,V0,THCK)
else:
    raise ValueError("Unrecognized potential type: %s" % POTENTIAL)

Vmax = V.max()    
dt   = hbar/(2*hbar**2/(m*dx**2)+Vmax)         
c1   = hbar*dt/(m*dx**2)                       
c2   = 2*dt/hbar                               
c2V  = c2*V  
print 'One-dimensional Schrodinger equation - time evolution'
print 'Wavepacket energy:   ',E
print 'Potential type:      ',POTENTIAL

print 'Barrier thickness:   ',THCK
psi_r = np.zeros((3,N)) #  Real
psi_i = np.zeros((3,N)) #  Imaginario
psi_p = np.zeros(N,)   #magnitud
                         
PA = 0                 
PR = 1                 
FU = 2                


xn = range(1,N/2)
x = X[xn]/dx    
gg = Gaussian(x,x0,sigma)
cx = np.cos(k0*x)
sx = np.sin(k0*x)
psi_r[PR,xn] = cx*gg
psi_i[PR,xn] = sx*gg
psi_r[PA,xn] = cx*gg
psi_i[PA,xn] = sx*gg
psi_p = psi_r[PR]**2 + psi_i[PR]**2

P   = dx * psi_p.sum()                     
nrm = np.sqrt(P)
psi_r /= nrm
psi_i /= nrm
psi_p /= P

pylab.figure()
xmin = X.min()
xmax = X.max()
ymax = 1.5*(psi_r[PR]).max()
pylab.axis([xmin,xmax,-ymax,ymax])

lineR, = pylab.plot(X,psi_r[PR],'b',alpha=0.7,label='Real')
lineI, = pylab.plot(X,psi_i[PR],'r',alpha=0.7,label='Imag')
lineP, = pylab.plot(X,6*psi_p,'k',label='Prob')
pylab.title('Grafica - Particula Libre', ha="center", fontsize=15,color="k", fontweight = 'bold')

if Vmax !=0 :
    s
    Efac = ymax/2.0/Vmax
    V_plot = V*Efac
    pylab.plot(X,V_plot,':k',zorder=0)   
    fillax(X,V_plot, facecolor='y', alpha=0.2,zorder=0)
   
    pylab.axhline(E*Efac,color='g',label='Energy',zorder=1)
pylab.legend(loc='lower right')
pylab.draw()

pylab.xlim(xmin,xmax)

IDX1 = range(1,N-1)                           
IDX2 = range(2,N)                             
IDX3 = range(0,N-2)                           
for t in range(T+1):
  
    psi_rPR = psi_r[PR]
    psi_iPR = psi_i[PR]

    psi_i[FU,IDX1] = psi_i[PA,IDX1] + \
                      c1*(psi_rPR[IDX2] - 2*psi_rPR[IDX1] +
                          psi_rPR[IDX3])
    psi_i[FU] -= c2V*psi_r[PR]

    psi_r[FU,IDX1] = psi_r[PA,IDX1] - \
                      c1*(psi_iPR[IDX2] - 2*psi_iPR[IDX1] +
                          psi_iPR[IDX3])
    psi_r[FU] += c2V*psi_i[PR]

    psi_r[PA] = psi_rPR
    psi_r[PR] = psi_r[FU]
    psi_i[PA] = psi_iPR
    psi_i[PR] = psi_i[FU]

    if t % Tp == 0:
        
        psi_p = psi_r[PR]**2 + psi_i[PR]**2
  
        lineR.set_ydata(psi_r[PR])
        lineI.set_ydata(psi_i[PR])
       
        lineP.set_ydata(15*psi_p)

        pylab.draw()

pylab.ioff()
pylab.show()
