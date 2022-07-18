'''
Created on May 16, 2022

@author: 31417
'''
# coding=utf8
from pylab import *
from numpy import *
from N2fixBudget import *
from E_BiosynthesisFromNH4 import *
from matplotlib.pyplot import *
from FigSetting2 import *
# this code is for the trichomes number supplementary figure

def QcFromV(V): #Compute Qc from V #Menden-Deuer 2000
    return 0.216*V**0.939/12 #(pmol C cell-1)

def QcFromVdiatom(V):
    return 0.288*V**0.811/12 #(pmol C cell-1)

def QnFromV(V): #Compute Qn from V #Menden-Deuer 2000
    return 0.118*V**0.849/14 #(pmol N cell-1)

# define a function use number of trichomes as the only variable, assume one heterocyst for one trichome
def comp(Nh):
    CN = 6.6
    
    Nv = 4 #Number of the vegetative cells per heterocysts
    #Nh = 2 #Number of the heterocysts per host diatom
    
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
    Mu = 0.51 #(d-1) Growth rate (Follett 2018)
    am = array([0,0.01,0.02,0.03,0.04,0.0467])# from no NH4+ uptake to all N from NH4+ uptake
    km = 0.483 # half saturate concentration (mmol/m-3)
    vmax = 1.16 # maximum uptake rate (pmol N C-1 d-1)
    Fn2fixN = Mu*(QnV + QnH + QnD)-(1.16*QcD)*(am/(am+0.483)) #N2 fixation rate (pmol N cell-1 d-1)
    vnh4 = (1.16*QcD)*(am/(am+0.483)) # NH4 uptake rate (pmol N cell-1 d-1)
    
    # N fixation can only happen when needed, if it is not needed, it should be 0 
    for i in range(size(Fn2fixN)):    
        if Fn2fixN[i] >= 0:
            Fn2fixN[i] = Fn2fixN[i] 
        else:
            Fn2fixN[i] = 0
        
    N_total = Fn2fixN+vnh4 # total N supply (pmol N cell-1 d-1) 
        
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
    return VD, DV_trans, FphoD, FphoV, Fn2fixN, vnh4, Fpho

N_t1,C_t1,FphoD01,FphoV01,NF1,NH1,Fpho1 = comp(1)
N_t2,C_t2,FphoD02,FphoV02,NF2,NH2,Fpho2 = comp(2)
N_t3,C_t3,FphoD03,FphoV03,NF3,NH3,Fpho3 = comp(3)
N_t4,C_t4,FphoD04,FphoV04,NF4,NH4,Fpho4 = comp(4)
N_t5,C_t5,FphoD05,FphoV05,NF5,NH5,Fpho5 = comp(5)

# NH4+ concentration gradient    
am = array([0,0.01,0.02,0.03,0.04,0.0467])

#'#FFCC00','#DF9900','#C06600','#A03300','#800000'
#thistle,plum,violet,m,purple

figure("Transfer N between Dia and Tricho")
plot(am,comp(1)[0],label = "trichomes*1",color = "#FFCC00")
plot(am,comp(2)[0],label = "trichomes*2",color = "#DF9900")
plot(am,comp(3)[0],label = "trichomes*3",color = "#C06600")
plot(am,comp(4)[0],label = "trichomes*4",color = "#A03300")
plot(am,comp(5)[0],label = "trichomes*5",color = "#800000")
plot(am,array([0,0,0,0,0,0]),linestyle = "dashed", color = "Black")
xlabel('NH$_{4}$$^{+}$ (mmol m$^{-3}$)')
ylabel('N Transfer (%)')
title("N transfer from Trichomes to Diatom")
legend()
savefig('figure2_1.png',dpi=300)

figure("Transfer C between Dia and Tricho")
plot(am,comp(1)[1],label = "trichomes*1",color = "#FFCC00" )
plot(am,comp(2)[1],label = "trichomes*2",color = "#DF9900")
plot(am,comp(3)[1],label = "trichomes*3",color = "#C06600")
plot(am,comp(4)[1],label = "trichomes*4",color = "#A03300")
plot(am,comp(5)[1],label = "trichomes*5",color = "#800000")
xlabel('NH$_{4}$$^{+}$ (mmol m$^{-3}$)')
ylabel('C Transfer (%)')
title('C transfer from Diatom to Trichomes')
legend()
savefig('figure2_2.png',dpi=300)


figure("C source from D")
plot(am,comp(1)[2],label = "trichomes*1",color = "#FFCC00" )
plot(am,comp(2)[2],label = "trichomes*2",color = "#DF9900")
plot(am,comp(3)[2],label = "trichomes*3",color = "#C06600")
plot(am,comp(4)[2],label = "trichomes*4",color = "#A03300")
plot(am,comp(5)[2],label = "trichomes*5",color = "#800000")
xlabel('NH$_{4}$$^{+}$ (mmol m$^{-3}$)')
ylabel('Photosynthesis rate (pmol C cell$^{-1}$ d$^{-1}$)')
title('Diatom Photosynthesis')
legend()
savefig('figure2_3.png',dpi=300)

figure("C source from V")
plot(am,comp(1)[3],label = "trichomes*1",color = "#FFCC00" )
plot(am,comp(2)[3],label = "trichomes*2",color = "#DF9900")
plot(am,comp(3)[3],label = "trichomes*3",color = "#C06600")
plot(am,comp(4)[3],label = "trichomes*4",color = "#A03300")
plot(am,comp(5)[3],label = "trichomes*5",color = "#800000")
xlabel('NH$_{4}$$^{+}$ (mmol m$^{-3}$)')
ylabel('Photosynthesis rate (pmol C cell$^{-1}$ d$^{-1}$)')
title('Trichomes Photosynthesis')
legend()
savefig('figure2_6.png',dpi=300)

figure("C source total")
plot(am,comp(1)[6],label = "trichomes*1",color = "#FFCC00" )
plot(am,comp(2)[6],label = "trichomes*2",color = "#DF9900")
plot(am,comp(3)[6],label = "trichomes*3",color = "#C06600")
plot(am,comp(4)[6],label = "trichomes*4",color = "#A03300")
plot(am,comp(5)[6],label = "trichomes*5",color = "#800000")
xlabel('NH$_{4}$$^{+}$ (mmol m$^{-3}$)')
ylabel('Photosynthesis rate (pmol C cell$^{-1}$ d$^{-1}$)')
title('Total Photosynthesis')
legend()
savefig('figure2_5.png',dpi=300)
#'#FFCC00','#DF9900','#C06600','#A03300','#800000'

figure("N2 fixation")
plot(am,comp(1)[4],label = "trichomes*1",color = "#FFCC00")
plot(am,comp(2)[4],label = "trichomes*2",color = '#DF9900')
plot(am,comp(3)[4],label = "trichomes*3",color = '#C06600')
plot(am,comp(4)[4],label = "trichomes*4",color = '#A03300')
plot(am,comp(5)[4],label = "trichomes*5",color = '#800000')
plot(am,comp(5)[5],label = "NH$_{4}$$^{+}$ uptake", color = "green")# also show NH4+ uptake in the same figure, NH4+ uptake is a singlue line because of the below figure
xlabel('NH$_{4}$$^{+}$ (mmol m$^{-3}$)')
ylabel('N rate (pmol N cell$^{-1}$ d$^{-1}$)')
title('N$_{2}$ fixation and NH$_{4}$$^{+}$ uptake')
legend()
savefig('figure2_4.png',dpi=300)

figure("NH4+ uptake")
plot(am,comp(1)[5],label = "1")
plot(am,comp(2)[5],label = "2")
plot(am,comp(3)[5],label = "3")
plot(am,comp(4)[5],label = "4")
plot(am,comp(5)[5],label = "5")
xlabel('NH4+ mmol m-3')
ylabel('N (pmol N cell-1 d-1)')
legend()


show()



