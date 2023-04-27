import numpy as np
import matplotlib.pyplot as plt
from math import exp, ceil, floor
from scipy import integrate
from scipy import interpolate
from scipy.optimize import curve_fit
from Sr327_loaddata import combineDATAnLARA, rhofactor
import sympy
import pandas as pd

if __name__ == "__main__":
    j_path = '../../Data/Sr327/Julian_Tom_data/'
    d_path = '../../Data/Sr327/Di_Tom_data/'
    cxcombined = combineDATAnLARA(j_path+'Version4Data/coolingdownP1',j_path+'Version4Data/coolingdownP1signal')
    cxcombinedtwo = combineDATAnLARA(j_path+'Version4Data/coolingdownP2',j_path+'Version4Data/coolingdownP2signal')
    cxcombinedfour = combineDATAnLARA(j_path+'Version4Data/heatingupP1',j_path+'Version4Data/heatingupP1signal')
    ndata = pd.concat([cxcombined[(cxcombined['Temperature']>146)&(cxcombined['Temperature']<292)]+0.0000015,cxcombinedfour[(cxcombinedfour['Temperature']>=7)&(cxcombinedfour['Temperature']<=146)].iloc[::-1],cxcombined[cxcombined['Temperature']<7],cxcombinedtwo])

    data = [[loadoldDATAlockin(d_path+'Sr327_Dec_05_2017_1_new.lvm'),'original c-axis',1.6],[loadDATAlockin(j_path+'Di_sample/isolation_c_cooling_down_2'),'c-axis',1.6],[loadDATAlockin(j_path+'Di_sample/isolation_d_coolingdown'),'diagonal',1.6],[loadDATAlockin(j_path+'Di_sample/isolation_x_coolingdown'),'in-plane',1.6],[ndata,'long c-axis',100*10**(-3)*rhofactor(5,1200,A,L)/83.35]]