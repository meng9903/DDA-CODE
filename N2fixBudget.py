'''
Created on 2013/01/05


#----------------------------------------------------------------------------------------------
From 617_02_01_ATP_50_reflected.py
#----------------------------------------------------------------------------------------------
@author: ag105020
'''
from pylab import *
from numpy import *
##########################################################################################################################################
def N2fixBudget(ep):   #ep: energy efficiency for the production of energy and the consumption of energy

    C=6             #number of C in one carbohydrate (ex. C(glucose)=6)
    n=6.6             #number of C in one bacterial biomass (BB)
    a=7             #number of H in one BB
    b=2             #number of O in one BB
    c=1             #number of N in one BB
    d=4*n+a-2*b-3*c #inverse number in the coefficient of BB in the half reaction for one e-
    BB=12*n+1*a+16*b+14*c #(g/mol): mass of BB (bacterial biomass)
    y=c/d           #coefficient of NH4+ in the half reaction of BB production
    z=1/4           #coefficient of NH4+ in the half reaction of nitrogen-fixation
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Computation of fe0, fpr and fn considering material, redox and energy balance
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #fs00=(1/y)/(1/y+1/z)        #The ratio of electron used for protein synthesis
    fs00 = 0
    #fn00=(1/z)/(1/y+1/z)        #The ration of electron used for nitrogen fixation
    fn00 = 1
    #ep=0.22                   #energy efficiency for the production of energy and the consumption of energy
    dgc0=41.35      #The free energy necessary for the half reaction of glucose production (kJ/e-mol)
    dgATP=50      #(kJ/ATP): energy produced by the reaction of ATP -> ADP (147-19)
    dgn=2*dgATP-dgc0*ep   #The free energy necessary (dg) for the half reaction of nitrogen fixation (kJ/e-mol)
    dgp=35.09-dgc0  #dg for production of pyruvate from glucose (kJ/e-mol))
    dgpc=3.33*1/d*(12*n+1*a+16*b+14*c)  #dg for the production of BB (bacterial biomass) from pyruvate) (147-17)
    dgr=-120.07     #-dg for the energy production pathway (kJ/e-mol)
    if dgp<0:       
        ep1=1/ep    #change ep1 depending on the sign of dgp    
    else:
        ep1=ep
    A=(fn00*dgn+fs00*(dgp/ep1+dgpc/ep))/(-ep*dgr) #A is related to fs0 and fe0
    fe0=A/(1+A)     #the ratio of electron used for energy production
    fs0=1-fe0       #the ratio of electron used for biomass synthesis+nitrogen fixation
    fpr=fs0*fs00    #the ratio of electron used for biomass synthesis
    fn=fs0*fn00     #the ratio of electron used for nitrogen fixation        
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #2.--Stoichiometry (to get E1(E for the case O2cri>O2in))------------------
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    S1=array([["","CH","H2O","CO2","O2","HCO3-","NH4+","N2","H2","BB","H+","e-"],
              ["'-Rd",1/24,0.25,-0.25,0.,0.,0.,0.,0.,0.,-1.,-1.],
              ["Ra",0.,-0.5,0.,0.25,0.,0.,0.,0.,0.,1.,1.],
              ["Rpr",0.,-(2*n-b+c)/d,(n-c)/d,0.,c/d,c/d,0.,0.,-1/d,1.,1.],
              ["Rn",0.,0.,0.,0.,0.,-0.25,0.125,-0.125,0.,1.25,1.]])
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S2 (f*R for electron acceptance)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Mu=array([[1],[fe0],[fpr],[fn]])  #column of f
    S2=copy(S1[1:5,1:12])           #use copy so S2 does not respond to the change in S1
    #S3=copy(S2[1:5,1:1,3])
    
    #convert S2 from string data type to float64 data type
    S2=S2.astype(float64)
    #add numbers for columns and raws for counting 
    S21=arange(1,12,1)
    S22=arange(0,5,1)
    S22=S22.reshape(5,1)
    S2=vstack((S21,S2))
    S2=hstack((S22,S2))
    #calculate f*R
    S2[1:5,1:12]=Mu*S2[1:5,1:12]
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S3 (f*R for electron donation)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    S3=vstack((S2[1,1:12],S2[1,1:12],S2[1,1:12]))
    S3=Mu[1:]*S3
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S4 (f*R for "electron donation + electron acceptance")
    # and RR, which is the entire reaction
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #print(S3)
    S4=S2[2:,1:]+S3     #S4 is f*R for "electron donation + electron acceptance"
    #print(S2[1:5,1:12])
    RR=S4[0]+S4[1]+S4[2]
    RR1=copy(RR)
    RR1=vstack((S1[0,1:],RR1))
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #output of each array into CSV files
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    savetxt("S1.csv", S1, delimiter=",",fmt='%s')
    savetxt("S2.csv", S2, delimiter=",",fmt='%2.3f')
    savetxt("S3.csv", S3, delimiter=",",fmt='%2.8f')
    savetxt("S4.csv", S4, delimiter=",",fmt='%2.8f')
    savetxt("RR.csv", RR, delimiter=",",fmt='%2.8f')
    savetxt("RR1.csv", RR1, delimiter=",",fmt='%s')

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #3.Calculation to obtain returning parameters
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    YcnN2fix = (S4[2,0]*C)/(S4[2,6]*2)
    YresN2fix = (S4[0,0]*C)/(S4[2,6]*2)
    
    return YcnN2fix,YresN2fix