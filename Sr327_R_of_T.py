import numpy as np
from math import exp



# c-axis resistivity, in mOhm-cm
# f(x) = 1.0/(exp((x-10)/10.) + 1)
# g(x) = 1.0/(exp((x-40)/20.) + 1)
# plot [0:300] (1 + 0.02*x**2)*f(x) + 0.08*x*g(x)*(1-f(x)) + (8.0+0.0005*x)*(1-g(x))*(1-f(x))

def R_c(T):  # in mu-Ohm-cm (hence factor of 1000.0)
  f = 1.0/(exp((T-10.0)/10.) + 1.)
  g = 1.0/(exp((T-40.0)/20.) + 1.)
  return 1000.0*(f*(1.0+0.02*T*T) + 0.08*T*g*(1.-f) + (8.0+0.0005*T)*(1.-g)*(1.-f))

# a-b plane resistivity, in mu-Ohm-cm
# f(x) = 1.0/(exp((x-20)/10.) + 1)
# plot [0:300] (1.7 + 0.03*x**2)*f(x) + 0.68*x*(1-f(x))

def R_ab(T): 
  f = 1.0/(exp((T-20.)/10.) + 1.0)
  return (1.7+0.03*T*T)*f + 0.68*T*(1.0-f)

