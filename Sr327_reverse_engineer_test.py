import numpy as np
from Outdated.Sr327_simulation import newR_c, R_ab
import matplotlib.pyplot as plt

def Mix(T,c,a):
    return 1/(1/(a*np.vectorize(R_ab)(T))+1/(c*np.vectorize(newR_c)(T)))
def Mix2(T,r):
    return ((1-r)*np.vectorize(R_ab)(T))+(r*np.vectorize(newR_c)(T))

T = np.linspace(0,300,100)
plt.plot(T,Mix(T,1,1500),color='blue')
plt.plot(T,Mix2(T,4/5),color='red')
# plt.plot(T,np.vectorize(R_ab)(T))
plt.plot(T,np.vectorize(newR_c)(T),color='black')
# plt.plot(T,np.vectorize(newR_c)(T)/np.vectorize(R_ab)(T))
plt.show()