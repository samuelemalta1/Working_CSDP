# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 17:25:00 2023

@author: Nick de Jong
"""
#Libraries
import numpy as np
import math
from scipy.optimize import fsolve
from scipy.optimize import least_squares
from scipy.optimize import minimize
import rocketcea as rc
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer
import Nozzle_turbine as NT
import Aux_classes as aux

#aux function, which calls the function for the selected cycle type
def TurboM(Default : aux.Default, prop : aux.Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
    # Default : [Default class object] contains information on the default values employed in the application
    # prop : [Propellant class object] contains the propellant data
    # p_a : [Pa] atmospheric pressure
    # Tf_cool : [K] [escalar] fuel temperature after cooling channel
    # dptcool : [Pa] [escalar] pressure losses in cooling channel
    # m : [kg/s] mass flow in the nozzle
    match Default.cycle_type:
        case 0: #Expander cycle
            turbo = EX(Default, prop, O_F, Tf_cool, dptcool, m)
        case 1: #Coolant bleed cycle
            turbo = CB(Default, prop, O_F, p_a, Tf_cool, dptcool, m)
        case 2: #Gas Generator
            turbo = GG(Default, prop, O_F, p_a, Tf_cool, dptcool, m)
        case 3: #Staged combustion cycle
            turbo = SC(Default, prop, O_F, p_a, Tf_cool, dptcool, m)
        case 4: #Electric motor to dirve pumps
            turbo = EL(Default, prop, O_F, Tf_cool, dptcool, m)
        case 5: #No turbomachinery
            dptvalve = Default.v_loss
            dptlines = 0.0
            dptmixing = 0.0
            return [((Default.p_to-dptvalve-dptlines)*O_F+Default.ptf-dptvalve-dptlines)/(1.0+O_F) - dptmixing, 0.0, 0.0, 0.0, 0]
        case 6: #Combustion tap off cycle
            turbo = TO(Default, prop, O_F, p_a, Tf_cool, dptcool, m) #Will not be currently implmented
        case _:
            print("Cycle Not recognized")
            return [0, 0, 0, 0, 1] #Returns 1 in last position, indicating aux to break
    
    turbo.results();
    return [turbo.ptinj, turbo.Wop, turbo.Wfp, turbo.Wt, turbo.br] #Returns 0 in last position, indicating aux to continue


#Function that computes the expander cycle
class EX:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_t : float #[-] efficiecny of turbine
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    pt1 : float #[Pa] Total pressure at turbine inlet
    pt2 : float #[Pa] Total pressure at turbine outlet
    ptinj : float #[Pa] Total pressure at injector inlet
    Wt : float #[W] Turbine power
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines

    #Flags
    br = False #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, Tf_cool : float, dptcool : float, m : float):
        print("Expander cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.dptvalve = DF.v_loss

    #Obtain results, calling optimization procedure and then computing variables of interest
    def results(self):
        #Minimization
        res = minimize(self.opt, [1.0e6], method = 'Nelder-Mead', bounds=[[0.0, 10.0e12]])
        self.dptfp = res["x"][0]

        #Computation of results
        self.dptop, self.pt1, self.pt2, self.ptinj = fsolve(self.equations,[1.0e6,1.0e7,1.0e6,1.0e6],self.dptfp)
        self.Wop = self.O_F/(self.O_F+1.0) * self.m * self.dptop / (self.eff_po*self.prop.o_dens)
        self.Wfp = 1.0/(self.O_F+1.0) * self.m * self.dptfp / (self.eff_pf*self.prop.f_dens_l)
        self.Wt = 1.0/(self.O_F+1.0) * self.m * self.eff_t * self.prop.fcp * self.Tf_cool * (1.0-(self.pt2/self.pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
        
        #Check to see if results are coherent
        if(abs(self.Wop+self.Wfp-self.Wt/self.eff_m) > 0.01 or self.Wop < 0.0 or self.Wfp < 0.0 or (not res["success"])):
            br = True

        #Debugging information
        print([self.dptop, self.dptfp, self.pt1, self.pt2, self.ptinj])
        print([self.Wop,self.Wfp,self.Wt])

    #Optimize for maximum chamber pressure
    def opt(self,dptfp):
        root = least_squares(self.equations,[1.0e6,1.0e7,1.0e6,1.0e6], args = dptfp, bounds = ((0.0,0.0,0.0,0.0),(10.0e10,10.0e10,10.0e10,10.0e10)))
        return -1.0*root["x"][3]

    #System of equations to be solved
    def equations(self,vars,dptfp):
        dptop, p1t, p2t, pinj = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - pinj,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - p1t,
            p2t - self.dptlines - pinj,
            self.O_F*dptop/(self.eff_po*self.prop.o_dens) + dptfp/(self.eff_pf*self.prop.f_dens_l) - self.eff_m*self.eff_t*self.prop.fcp*self.Tf_cool*(1.0-(p2t/p1t)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
            ]

    
#Function that computes staged combustion
class SC:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_t : float #[-] efficiecny of turbine
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling
    pa : float #[Pa] Ambient pressure

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    pt1 : float #[Pa] Total pressure at turbine inlet
    pt2 : float #[Pa] Total pressure at turbine outlet
    ptinj : float #[Pa] Total pressure at injector inlet
    Wt : float #[W] Turbine power
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power
    l : float #[-] fraction of oxidizer for pre-combustor
    T1t : float #[K] Temperature after pre-combustor

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines
    dptmix = 0.0 #[Pa] mixing total pressure losses
    dptcomb = 0.0 #[Pa] pressure losses in combustion

    #Flags
    br = False #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Staged Combustion cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.dptvalve = DF.v_loss

    #Obtain results, calling optimization procedure and then computing variables of interest
    def results(self):
        #Minimization
        res = minimize(self.opt, [1.0e5,0.01], method = 'Nelder-Mead', bounds=[[self.pa*1.2,10.0e12],[1.0e-5,1.0]])
        self.pt2 = res["x"][0]
        self.l = res["x"][1]

        #Computation of results
        self.dptop, self.pt1, self.dptfp, ptinj = fsolve(self.equations,[1.0e6,1.0e7,1.0e7,1.0e7],(self.pt2,self.l))
        self.Wop = self.O_F/(self.O_F+1.0) * self.m * self.dptop / (self.eff_po*self.prop.o_dens)
        self.Wfp = 1.0/(self.O_F+1.0) * self.m * self.dptfp / (self.eff_pf*self.prop.f_dens_l)
        self.Wt = (1.0+self.l*self.O_F)/(self.O_F+1.0) * self.m * self.eff_t * self.prop.fcp * self.T1t * (1.0-(self.pt2/self.pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
        
        #Check to see if results are coherent
        if(abs(self.Wop+self.Wfp-self.Wt/self.eff_m) > 0.01 or self.Wop < 0.0 or self.Wfp < 0.0 or (not res["success"])):
            br = True

        #Debugging information
        print([self.dptop, self.dptfp, self.pt1, self.pt2, self.ptinj])
        print([self.Wop,self.Wfp,self.Wt])

    #Optimize for maximum chamber pressure
    def opt(self,vars):
        root = least_squares(self.equations,[1.0e6,1.0e7,1.0e7,1.0e7], args = vars, bounds = ((10.0,self.pa,10.0,self.pa*1.5),(10.0e10,10.0e10,10.0e10,10.0e10)))
        Isp_m = self.get_Isp_m(root["x"][1]) #ver como pasar aqui la composición modificada.....
        return -Isp_m
    
    def get_Isp_m(self,pinj): #temporaty, needs modification
        return NT.Turbine_nozzle(self.m,pinj,self.prop,self.pa,self.df,self.prop.h_fuel,self.prop.h_ox,self.prop.f_dens_g,self.prop.o_dens)

    #System of equations to be solved
    def equations(self,vars,p2t,l):
        dptop, p1t, dptfp, pinj = vars
        pcomb = self.ptankf - self.dptvalve + dptfp - self.dptcool
        self.T1t = 1.0 #Function to get combustion temperature in pre-combustor
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - pinj,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - self.dptmix - self.dptcomb - p1t,
            p2t - self.dptlines - pinj,
            self.O_F*dptop/(self.eff_po*self.prop.o_dens) + dptfp/(self.eff_pf*self.prop.f_dens_l) - (1.0+l*self.O_F)*self.eff_m*self.eff_t*self.prop.fcp*self.T1t*(1.0-(p2t/p1t)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
            ]


#Function that computes coolant bleed
class CB:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_t : float #[-] efficiecny of turbine
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling
    pa : float #[Pa] Ambient pressure

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    pt1 : float #[Pa] Total pressure at turbine inlet
    pt2 : float #[Pa] Total pressure at turbine outlet
    ptinj : float #[Pa] Total pressure at injector inlet
    Wt : float #[W] Turbine power
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power
    l : float #[-] fraction of fuel for bleed

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines
    m_O = 1.0 #[kg/s] Oxidizer mass flow
    m_F = 1.0 #[kg/s] Fuel mass flow

    #Flags
    br = False #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Coolant bleed cycle selected")
        self.df = DF
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.m_O = O_F/(O_F+1.0) * m
        self.m_F = 1.0/(O_F+1.0) * m
        self.dptvalve = DF.v_loss

    #Obtain results, calling optimization procedure and then computing variables of interest
    def results(self):
        #Minimization
        res = minimize(self.opt, [1.0e5,0.01], method = 'Nelder-Mead', bounds=[[self.pa*1.2,10.0e12],[1.0e-5,1.0]])
        self.pt2 = res["x"][0]
        self.l = res["x"][1]

        #Computation of results
        self.dptop, self.pt1, self.dptfp = fsolve(self.equations,[1.0e6,1.0e7,1.0e7],(self.pt2,self.l))
        self.ptinj = self.pt1
        self.Wop = self.m_O * self.dptop / (self.eff_po*self.prop.o_dens)
        self.Wfp = (1.0/(1.0-self.l)) * self.m_F * self.dptfp / (self.eff_pf*self.prop.f_dens_l)
        self.Wt = (self.l/(1.0-self.l)) * self.m_F * self.eff_t * self.prop.fcp * self.Tf_cool * (1.0-(self.pt2/self.pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
        
        #Check to see if results are coherent
        if(abs(self.Wop+self.Wfp-self.Wt/self.eff_m) > 0.01 or self.Wop < 0.0 or self.Wfp < 0.0 or (not res["success"])):
            br = True

        #Debugging information
        print([self.dptop, self.dptfp, self.pt1, self.pt2, self.l])
        print([self.Wop,self.Wfp,self.Wt])

    #Optimize for maximum chamber pressure
    def opt(self,vars):
        root = least_squares(self.equations,[1.0e6,1.0e7,1.0e7], args = vars, bounds = ((10.0,self.pa,10.0),(10.0e10,10.0e10,10.0e10)))
        Isp_m = self.get_Isp_m(root["x"][1])
        Isp_a = self.get_Isp_a(root["x"][1], vars[0], vars[1])
        return -(Isp_m + Isp_a*(vars[1]/(1.0-vars[1])) * 1.0/(self.O_F+1.0))/(1.0+(vars[1]/(1.0-vars[1])) * 1.0/(self.O_F+1.0))
    
    #Get Isp of aux nozzle
    def get_Isp_m(self,pinj): #temporaty, needs modification
        return NT.Turbine_nozzle(self.m,pinj,self.prop,self.pa,self.df,self.prop.h_fuel,self.prop.h_ox,self.prop.f_dens_g,self.prop.o_dens,self.O_F)
    
    #get Isp of open cycle auxiliary nozzle
    def get_Isp_a(self,pt1,pt2,l):
        T2t = self.Tf_cool*(pt2/pt1)**((self.prop.f_gamma-1.0)/self.prop.f_gamma)
        return (l/(1.0-l))*self.m_F*math.sqrt(2.0*self.prop.R_f*T2t*self.prop.f_gamma/(self.prop.f_gamma-1.0) * (1.0 - (self.pa/pt2)**((self.prop.f_gamma-1.0)/self.prop.f_gamma)))

    #System of equations to be solved
    def equations(self,vars,p2t,l):
        dptop, p1t, dptfp = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - p1t,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - p1t,
            self.O_F*dptop/(self.eff_po*self.prop.o_dens) + (1.0/(1.0-l))*dptfp/(self.eff_pf*self.prop.f_dens_l) - (l/(1.0-l))*self.eff_m*self.eff_t*self.prop.fcp*self.Tf_cool*(1.0-(p2t/p1t)**((self.prop.f_gamma-1.0)/self.prop.f_gamma))
            ]


#Function that computes Gas Generator
class GG:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_t : float #[-] efficiecny of turbine
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling
    pa : float #[Pa] Ambient pressure

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    pt1 : float #[Pa] Total pressure at turbine inlet
    pt2 : float #[Pa] Total pressure at turbine outlet
    ptinj : float #[Pa] Total pressure at injector inlet
    Wt : float #[W] Turbine power
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines
    m_O = 1.0 #[kg/s] Oxidizer mass flow
    m_F = 1.0 #[kg/s] Fuel mass flow

    #Flags
    br = False #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Gas generator cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.dptvalve = DF.v_loss


#Function that computes Combustion Tap-Off cycle
class TO:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_t : float #[-] efficiecny of turbine
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling
    pa : float #[Pa] Ambient pressure

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    pt1 : float #[Pa] Total pressure at turbine inlet
    pt2 : float #[Pa] Total pressure at turbine outlet
    ptinj : float #[Pa] Total pressure at injector inlet
    Wt : float #[W] Turbine power
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power
    l : float #[-] fraction of fuel for bleed

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines
    m_O = 1.0 #[kg/s] Oxidizer mass flow
    m_F = 1.0 #[kg/s] Fuel mass flow

    #Flags
    br = False #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, p_a : float, Tf_cool : float, dptcool : float, m : float):
        print("Tap-off cycle not supported")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_t = DF.Eff_t
        self.eff_m = DF.Eff_m
        self.pa = p_a
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.dptvalve = DF.v_loss


#Function that computes Electric driven pumps cycle
class EL:
    #Global input
    ptanko : float #[Pa] Oxidizer tank pressure
    ptankf : float #[Pa] Fuel tank pressure
    prop : aux.Propellant #[Object from Propellant class] contains propellant data
    O_F : float #[-] Oxidizer to fuel ratio
    eff_po : float #[-] efficiency of oxidizer pump
    eff_pf : float #[-] efficiency of fuel pump
    eff_m : float #[-] mechanical efficiency of turbine-pump coupling
    wmotor : float #[W] Power of electrical motor
    Wt : float #[W] Equal to wmotor, used for conveniance

    #Software input
    Tf_cool : float #[K] Temperature of fuel after cooling channel
    dptcool : float #[Pa] Total pressure losses in cooling channel
    m : float #[kg/s] mass flow in nozzle

    #Output
    dptop : float #[Pa] Total pressure increase in oxidizer pump
    dptfp : float #[Pa] Total pressure increase in fuel pump
    ptinj : float #[Pa] Total pressure at injector inlet
    Wop : float #[W] Oxidizer pump power
    Wfp : float #[W] Fuel pump power

    #Auxiliary
    dptvalve = 0.0 #[Pa] Total pressure losses in valves
    dptlines = 0.0 #[Pa] Total pressure losses in lines

    #Flags
    br = False #If true the aux program breaks after this function (an error has occured)

    #Initialize values
    def __init__(self, DF : aux.Default, prop : aux.Propellant, O_F : float, Tf_cool : float, dptcool : float, m : float):
        print("Coolant bleed cycle selected")
        self.ptanko = DF.p_to
        self.ptankf = DF.ptf
        self.prop = prop
        self.O_F = O_F
        self.eff_po = DF.Eff_po
        self.eff_pf = DF.Eff_pf
        self.eff_m = DF.Eff_m
        self.wmotor = DF.Wmotor
        self.Wt = DF.Wmotor
        self.Tf_cool = Tf_cool
        self.dptcool = dptcool
        self.m = m
        self.dptvalve = DF.v_loss

    #Obtain results, no optimization procedure needed for this cycle
    def results(self):
        #Computation of results
        self.dptop, self.dptfp, self.ptinj = fsolve(self.equations,[1.0e6,1.0e6,1.0e6])
        self.Wop = self.O_F/(self.O_F+1.0) * self.m * self.dptop / (self.eff_po*self.prop.o_dens)
        self.Wfp = 1.0/(self.O_F+1.0) * self.m * self.dptfp / (self.eff_pf*self.prop.f_dens_l)

        #Check to see if results are coherent
        if(abs(self.Wop+self.Wfp-self.Wt/self.eff_m) > 0.01 or self.Wop < 0.0 or self.Wfp < 0.0):
            br = True

        #Debugging information
        print([self.dptop, self.dptfp, self.ptinj])
        print([self.Wop,self.Wfp,self.wmotor])

    #Equations to be solved
    def equations(self,vars):
        dptop, dptfp, pinj = vars
        return [
            self.ptanko - self.dptvalve + dptop - self.dptlines - pinj,
            self.ptankf - self.dptvalve + dptfp - self.dptcool - pinj,
            self.O_F*dptop/(self.eff_po*self.prop.o_dens) + dptfp/(self.eff_pf*self.prop.f_dens_l) - self.eff_m*self.wmotor
            ]


#Values used during tests, corresponding to Global or software input
O_F_ = 5.0 #[-]
Pa_ = 1.0e5 #[Pa] ambient pressure
Tf_cool_ = 500.0 #[K] temperature after cooling
dptcool_ = 1.0e5 #[Pa] pressure drop in cooling channel
m_ = 2.0 #[kg/s] mass flow in the nozzle
default = aux.Default(0)
prop = aux.Propellant(0)

#main Function
if __name__ == '__main__':
    print('Loading...')
    default.cycle_type = 0
    print(TurboM(default, prop, O_F_, Pa_, Tf_cool_, dptcool_, m_))
    print('\nProcess Terminated')