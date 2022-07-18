'''
Created on May 16, 2022

@author: 31417
'''
from pylab import *
from numpy import *
from N2fixBudget import *
from E_BiosynthesisFromNH4 import *
from matplotlib.pyplot import *
import warnings
warnings.filterwarnings("ignore")
from FigSetting2 import *

# The code is to calculate the values in the flux figure
def QcFromV(V): #Compute Qc from V #Menden-Deuer 2000
    return 0.216*V**0.939/12 #(pmol C cell-1)

def QcFromVdiatom(V):
    return 0.288*V**0.811/12 #(pmol C cell-1)

def QnFromV(V): #Compute Qn from V #Menden-Deuer 2000
    return 0.118*V**0.849/14 #(pmol N cell-1)
# define a function with growth rate and NH4+ concentration as input data
def comp(Mu,am):
    CN = 6.6
    
    Nv = 4 #Number of the vegetative cells per heterocysts
    Nh = 2 #Number of the heterocysts per host diatom
    
    R = 8
    #R = arange(6,11+0.1,0.1)
    
    Vv0 = 18.845 #(um3) Cell Volume of the vegetative cells from "03 Cell volume.xlsx"
    Vh0 = 61.0125 #(um3) Cell volume of the heterocysts "03 Cell volume.xlsx"
    Vd00 = 3493.529 #(um3) Cell volume of the host diatom "03 Cell volume.xlsx"
    R0 = (3*Vd00/4/pi)**(1/3)
    R0x = (R0,R0)
    R0y = (0,1400)
    Vd0 = 4/3*pi*R**3
    
    Vv = Vv0*Nv*Nh #(um3) Cell Volume of the vegetative cells from "03 Cell volume.xlsx"
    Vh = Vh0*Nh #(um3) Cell volume of the heterocysts "03 Cell volume.xlsx"
    Vd = Vd0 #(um3) Cell volume of the host diatom "03 Cell volume.xlsx"
    
    QcV0 = QcFromV(Vv0) #(pmol C cell-1)
    QcH0 = QcFromV(Vh0) #(pmol C cell-1)
    QcD0 = QcFromVdiatom(Vd0) #(pmol C cell-1)
    
    QcV = QcV0*Nv*Nh
    QcH = QcH0*Nh
    QcD = QcD0
    
    QcC = QcV + QcH #C as "Chain"
    Qc = QcV + QcH + QcD
    
    QnV = QcV/CN #(pmol N cell-1)
    QnH = QcH/CN #(pmol N cell-1)
    QnD = QcD/CN #(pmol N cell-1)
    
    QnC = QnV + QnH
    Qn = QnV + QnH + QnD
    
 
    Fn2fixN = Mu*(QnV + QnH + QnD)-(1.16*QcD)*(am/(am+0.483)) #N fixation rate (pmol N cell-1 d-1)
    vnh4 = (1.16*QcD)*(am/(am+0.483)) # NH4+ uptake rate (pmol N cell-1 d-1)
    
    # N fixation can only happen when needed, if it is not needed, it should be 0
    if Fn2fixN >= 0:
            Fn2fixN = Fn2fixN 
    else:
            Fn2fixN = 0
        
    N_total = Fn2fixN+vnh4 # total N supple (pmol N cell-1 d-1)
        
    E = E_BiosynthesisFromNH4(0.6) #(mol C mol C-1) CO2 production rate
    
    YcnN2fix,YresN2fix = N2fixBudget(0.6) #(mol C mol N-1) Converting N2 fixation to respiration and C consumption for electron
    
    Fn2fixC = Fn2fixN*YcnN2fix + Fn2fixN*YresN2fix #(pmol C cell-1 d-1) C consumption for N2 fixation
    
    Fpho = Mu*(QcV + QcH + QcD)*(1+E) + Fn2fixC #(pmol C cell-1 d-1) Photosynthesis rate
    
    FphoV = Fpho*QnV/(QnV + QnD) #C percentage from vegetative cells photosynthesis
    FphoD = Fpho*QnD/(QnV + QnD) #C percentage from diatom photosynthesis
    D_growth = QcD*Mu/Fpho*100 # C percentage used in diatom growth
    D_res = QcD*Mu*E/Fpho*100 # C percentage used in diatom respiration
    DV_trans = (FphoD-QcD*Mu-QcD*Mu*E)/Fpho*100 #C percentage used in C transfer
    
    Nc_fix = Fn2fixC/Fpho*100 # C percentage used in N2 fixation 

    Nfix_p = Fn2fixN/N_total*100 # N percentage from N2 fixation
    NH4_p = vnh4/N_total*100 # N percentage from NH4+ uptake
    DN_growth = Mu*QnD/N_total*100 # N percentage used in diatom growth
    VN_growth = Mu*QnV/N_total*100 # N percentage used in vegatitive growth
    HN_growth = Mu*QnH/N_total*100 # N percentage used in heterocysts growth
    d = HN_growth+VN_growth # N percentage used in trichomes biosynthesis
    
    tot = Mu*QcD*(1+E)+Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC # all C used
    a = Mu*QcD*(1+E)/tot*100 # C portion used in diatom biosynthesis
    b = (Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E)/tot*100 # C portion used in trichomes biosynthesis
    c = Fn2fixC/tot*100 # C portion used in N fixation
    
    # Calculate N transfer
    for i in range(size(am)):
        if Nfix_p> (QnV+QnH)/(QnD+QnV+QnH):
            VD = Nfix_p - HN_growth - VN_growth
        else:
            VD = -(NH4_p-DN_growth)
        
    # Calculate N transfer NH4 and the NH4 supply for all     
    am_0 = ((Mu*QnD*0.483)/(1.16*QcD)/(1-(Mu*QnD)/(1.16*QcD)))
    am_1 = (0.483*Mu*(QnV + QnH + QnD))/(1.16*QcD-Mu*(QnV + QnH + QnD)) 
    return FphoD/Fpho*100, FphoV/Fpho*100,a,DV_trans,b,c,Nfix_p,NH4_p,VD,HN_growth+VN_growth,DN_growth


# print flux values (unit %)
print(comp(0.4,0.01))
print(comp(0.8,0.01))
print(comp(0.4,0.036))
print(comp(0.8,0.036))


