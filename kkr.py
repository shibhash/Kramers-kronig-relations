import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import scipy.constants as sc


wp=1.35*10e15
w0=0.6*10e15
w = np.arange(0, 1.4*10e15, 0.01*10e15) #angular frequency range 
gamma=0.15*10e15

def imaginary_e(w):     #complex valued function obtained for the dielectric permittivity in the Lorentz model
    x=w0**2
    y=w**2
    z=x-y
    e_i= (wp**2)*gamma*w/(z**2+w**2*gamma**2)
    return e_i;
    
def real_e(w):         #Real valued function obtained for the dielectric permittivity in the Lorentz model
    x=w0**2
    y=w**2
    z=x-y
    e_r= 1.0 + (wp**2)*z/(z**2+w**2*gamma**2)
    return e_r;
    
eps=sc.epsilon_0


def integration(w):     # Principal value calcualation with delta=100
    integrand = lambda w1 : w1*(imaginary_e(w1))/((w1**2-w**2))
    integral1, err1 = quad(integrand, w+100, 5*10e15) 
    integral2, err2 = quad(integrand, 0, w-100)
    return eps*(1+(2/(np.pi*eps))*(integral1+integral2));


KK_real=[]              #assigning values to array
for i in range(len(w)):
    KK_real.append(integration(w[i]))
            


#plt.plot(w,imaginary_e(w), linestyle='--') 
line1, =plt.plot(w, real_e(w), label='Lorentz model', linestyle='--') #plot
line2, =plt.plot(w, KK_real, label='Kramers-Kronig Relations', linewidth=3)
plt.xlabel('Angular frequency($\omega$)')
plt.ylabel('Re($\epsilon$)')
plt.title('Re($\epsilon$) vs. Angular frequency($\omega$)')
plt.grid(True)
first_legend = plt.legend(handles=[line1], loc=1)
ax = plt.gca().add_artist(first_legend)
plt.legend(handles=[line2], loc=4)
plt.show()
