import matplotlib.pyplot as plt
import numpy as np

el = 47*10**(3)
debye = 57
lambdap = 18.5

print(debye)
print(lambdap)

def resistivity(T, lambdap, debye):
    if (T < debye/np.sqrt(lambdap)).all():
        ratio = np.exp(-lambdap/2)*(1+((np.pi**2 * lambdap)/3) *(T/debye)**2)
    elif ((debye/np.sqrt(lambdap) <= T)).all() and ((T < debye)).all():
        ratio = np.exp(-lambdap/2 + (lambdap/3)*((np.pi * T)/debye)**2)
    elif (( debye <= T )).all() and (( T < lambdap * debye )).all():
        ratio = np.sqrt((3* (np.pi ** 3) * T)/(4*lambdap * debye))*np.exp(-(lambdap*debye)/(12*T))
    elif (lambdap * debye <= T).all():
        ratio = 1- (lambdap/9)*(debye/T)
    else:
        print("ERROR!")
        ratio = 1
    return (el*ratio)**(-1)

T1 = np.linspace(0, debye/np.sqrt(lambdap),100)
T2 = np.linspace(debye/np.sqrt(lambdap), debye,100)
T3 = np.linspace(debye, lambdap * debye,100)
T4 = np.linspace(lambdap * debye, 4*lambdap * debye,100)

plt.plot(T1[:-1],resistivity(T1[:-1],lambdap,debye),label = "Region 1")
plt.plot(T2[:-1],resistivity(T2[:-1],lambdap,debye),label = "Region 2")
plt.plot(T3[:-1],resistivity(T3[:-1],lambdap,debye),label = "Region 3")
plt.plot(T4,resistivity(T4,lambdap,debye),label = "Region 4")
plt.legend()
plt.xlim(2,300)
plt.ylim(-1,100)
plt.show()