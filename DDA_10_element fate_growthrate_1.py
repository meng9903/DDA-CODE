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

# this code is for the growth rate figures both in main text and in supplementary

def QcFromV(V): #Compute Qc from V #Menden-Deuer 2000
    return 0.216*V**0.939/12 #(pmol C cell-1)

def QcFromVdiatom(V):
    return 0.288*V**0.811/12 #(pmol C cell-1)

def QnFromV(V): #Compute Qn from V #Menden-Deuer 2000
    return 0.118*V**0.849/14 #(pmol N cell-1)

# define a function use growth rate as the only variable, assume one heterocyst for one trichome
def comp(Mu):
    CN = 6.6 # redfield ratio C:N
    
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
 
    
    #Maybe I can make it growth rate array and use it for the x-axis.
    
    am = arange(0,0.0467,0.001)# an array for NH4+ gradient
    Fn2fixN = Mu*(QnV + QnH + QnD)-(1.16*QcD)*(am/(am+0.483)) #N fixation rate (pmol N cell-1 d-1)
    vnh4 = (1.16*QcD)*(am/(am+0.483)) # NH4+ uptake rate (pmol N cell-1 d-1)
    
    # N fixation can only happen when needed, if it is not needed, it should be 0
    for i in range(size(Fn2fixN)):    
        if Fn2fixN[i] >= 0:
            Fn2fixN[i] = Fn2fixN[i] 
        else:
            Fn2fixN[i] = 0
        
    N_total = Fn2fixN+vnh4 # total N supple (pmol N cell-1 d-1)
        
    E = E_BiosynthesisFromNH4(0.6) #(mol C mol C-1) CO2 production rate
    
    YcnN2fix,YresN2fix = N2fixBudget(0.6) #(mol C mol N-1) Converting N2 fixation to respiration and C consumption for electron
    
    Fn2fixC = Fn2fixN*YcnN2fix + Fn2fixN*YresN2fix #(pmol C cell-1 d-1) C consumption for N2 fixation
    
    Fpho = Mu*(QcV + QcH + QcD)*(1+E) + Fn2fixC #(pmol C cell-1 d-1) Photosynthesis rate
    
    FphoV = Fpho*QnV/(QnV + QnD) #(pmol C cell-1 d-1) Photosynthesis rate by vegetative cells
    FphoD = Fpho*QnD/(QnV + QnD) #(pmol C cell-1 d-1) Photosynthesis rate by diatom cells
    DV_trans = (FphoD-QcD*Mu-QcD*Mu*E)/Fpho*100 # C percentage used in transfer from diatoms to vegetative cells

    Nfix_p = Fn2fixN/N_total*100 # N percentage from N2 fixation
    NH4_p = vnh4/N_total*100 # N percentage from NH4 uptake
    DN_growth = Mu*QnD/N_total*100 # N percentage used in diatom growth
    VN_growth = Mu*QnV/N_total*100 # N percentage used in vegetative cells growth
    HN_growth = Mu*QnH/N_total*100 # N percentage used in heterocysts growth
    VD = zeros(size(am))*nan # create a array for N transfer
    # Calculate N transfer
    for i in range(size(am)):
        if Nfix_p[i]> (QnV+QnH)/(QnD+QnV+QnH):
            VD[i] = Nfix_p[i] - HN_growth[i] - VN_growth[i]
        else:
            VD[i] = -(NH4_p[i]-DN_growth[i])
        
    # Calculate N transfer NH4 and the NH4 supply for all     
    am_0 = ((Mu*QnD*0.483)/(1.16*QcD)/(1-(Mu*QnD)/(1.16*QcD)))
    am_1 = (0.483*Mu*(QnV + QnH + QnD))/(1.16*QcD-Mu*(QnV + QnH + QnD))
    # Return the value we need 
    return VD, DV_trans, am_0, FphoD, FphoV, Fn2fixN, Fpho

M = arange(0.3,0.8,0.01)# an array for different growth rate
am = arange(0,0.0467,0.001)# an array for NH4+ concentration 

# Create arrays for the value we need
VD = zeros([size(M),size(am)])
DV_trans = zeros([size(M),size(am)])
am_0 = zeros([size(M),size(am)])
FphoD = zeros([size(M),size(am)])
FphoV = zeros([size(M),size(am)])
Fn2fixN = zeros([size(M),size(am)])
Fpho = zeros([size(M),size(am)])


for i in range(size(M)):
    VD[i,] = comp(M[i])[0]

for i in range(size(M)):
    DV_trans[i,] = comp(M[i])[1]
    
for i in range(size(M)):
    am_0[i,] = comp(M[i])[2] 

for i in range(size(M)):
    FphoD[i,] = comp(M[i])[3]        

for i in range(size(M)):
    FphoV[i,] = comp(M[i])[4]  

for i in range(size(M)):
    Fn2fixN[i,] = comp(M[i])[5]  
    
for i in range(size(M)):
    Fpho[i,] = comp(M[i])[6]        
      
    
# Make figures    
figure(0)
pcolormesh(am,M,VD,shading = 'flat')#cmap = 'RdBu'
colorbar()
ylabel("growth rate (d$^{-1}$)")
xlabel("NH$_{4}$$^{+}$ (mmol m$^{-3}$)")
title("N transfer from Trichomes to Diatom (%)")
savefig('figure3_1.png',dpi=300)

figure(1)
pcolormesh(am,M,DV_trans,shading = 'flat')
colorbar()
ylabel("growth rate (d$^{-1}$)")
xlabel("NH$_{4}$$^{+}$ (mmol m$^{-3}$)")
title("C transfer from Diatom to Trichomes (%)")
savefig('figure3_2.png',dpi=300)

figure(2)
pcolormesh(am,M,am_0,shading = 'flat')
colorbar()
ylabel("growth rate (d$^{-1}$)")
xlabel("NH$_{4}$$^{+}$ (mmol m$^{-3}$)")
title("No transfer NH$_{4}$$^{+}$ (mmol m$^{-3}$)")
savefig('figure3_3.png',dpi=300)


figure(6)
pcolormesh(am,M,FphoD,shading = 'flat')
colorbar()
ylabel("growth rate (d$^{-1}$)")
xlabel("NH$_{4}$$^{+}$ (mmol m$^{-3}$)")
title("Diatom photosynthesis (pmol C cell$^{-1}$ d$^{-1}$)")
savefig('figure3_4.png',dpi=300)

figure(7)
pcolormesh(am,M,FphoV,shading = 'flat')
colorbar()
ylabel("growth rate (d$^{-1}$)")
xlabel("NH$_{4}$$^{+}$ (mmol m$^{-3}$)")
title("Trichomes photosynthesis (pmol C cell$^{-1}$ d$^{-1}$)")
savefig('figure3_5.png',dpi=300)

figure(8)
pcolormesh(am,M,Fn2fixN,shading = 'flat')
colorbar()
ylabel("growth rate (d$^{-1}$)")
xlabel("NH$_{4}$$^{+}$ (mmol m$^{-3}$)")
title("N$_{2}$ fixation rate (pmol N cell$^{-1}$ d$^{-1}$)")
savefig('figure3_6.png',dpi=300)

figure(9)
pcolormesh(am,M,Fpho,shading = 'flat')
colorbar()
ylabel("growth rate (d$^{-1}$)")
xlabel("NH$_{4}$$^{+}$ (mmol m$^{-3}$)")
title("Total photosynthesis (pmol C cell$^{-1}$ d$^{-1}$)")
savefig('figure3_7.png',dpi=300)




show()
    
