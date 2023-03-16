import math as mth

class Materials:
    def __init__(self, material, density, yieldstress_l, Emod, OpTemp_u, k, cost):
        self.material = material
        self.density = density
        self.yieldstress_l = yieldstress_l
        self.Emod = Emod
        self.OpTemp_u = OpTemp_u
        self.k = k
        self.cost = cost
        

Rhenium                 =       Materials('Rhenium',                            21000.0,  2300.0e6,   471.0e9, 2200.0,   48.0,   938.0)
Al_7075_T6              =       Materials('Aluminium 7075 T6',                  2810.0,   505.0e6,    72.0e9,  373.15,   130.0,  3.97)
Ti6Al4V                 =       Materials('Tungsten & Aluminium alloy',         4428.0,   110.0e6,    114.0e9, 693.15,   6.7,    23.8)
Haynes_188              =       Materials('Nickel Alloy',                       8980.0,   464.0e6,    244.0e9, 1423.0,   10.8,   30.6)
Inc_X_750               =       Materials('Inconel-750, Nickel-Chromium Alloy', 8280.0,   850.0e6,    222.0e9, 1090.0,   12.0,   20.4)
Inc_600                 =       Materials('Inconel-600, Nickel-Chromium Alloy', 8470.0,   290.0e6,    206.0e9, 982.0,    15.9,   20.2)
Inc_718                 =       Materials('Inconel-718, Nickel-Chromium Alloy', 8190.0,   1036.0e6,   211.0e9, 825.0,    11.4,   16.6)
Inc_A_286               =       Materials('A-286_Nickel-Chromium Alloy',        7920.0,   7920.0e6,   201.0e9, 700.0,    23.9,   5.43)
Columbium_c103          =       Materials('Niobium (Colombium) - cold rolled',  8600.0,   550.0e6,    130.0e9, 1255.0,   0.54,   225.0)
Copper_structural       =       Materials('Copper',                             8940.0,   40.0e6,     128.0e9, 573.0,    398.0,  5.60) 
D6AC_Steel              =       Materials('D6AC Steel Alloy',                   7870.0,   145.0e6,    200.0e9, 774.0,    52.0,   1.12)

#Coating Materials:
Copper                  =       Materials('Copper coating',                     8940.0,   0,          0,          573.0,    390.0,  6.4)
Narloy_Z                =       Materials('Copper Alloy Coating',               9130.0,   0,          0,          740.0,    290.0,  6.4)
GRCop_84                =       Materials('GRPCop_84 Copper alloy coating',     8756.0,   0,          0,          1000.0,   351.0,  6.4)
Silica                  =       Materials('Scilica Coating',                    1700.0,   0,          0,          1200.0,   0.55,   6.4)
Carbon                  =       Materials('Carbon-Carbon Matrix coating',       1950.0,   0,          0,          2400.0,   37.4,   6.4)


# ##Computing Mass nozzle: 
# #Nozzle Surface Fucntion:
# def Nozzle_surface(x,R):
#     total_dist = 0
#     for i in range(len(x)-1):
#         total_dist += mth.dist([x[i+1],R[i+1]],[x[i],R[i]])
#     vol = total_dist*2*mth.pi
#     return vol

class ReferenceEngine:
    def __init__(self, pc, thrust, arear, rt, mprop, rhoprop, FS, Material_NCG, Material_P, Material_V, mfrac_tube, mfrac_manifold, mfrac_jacket, mfrac_chamber, mfrac_pump, mfrac_valve, RefMass):
        self.pc = pc
        self.thrust = thrust
        self.arear = arear
        self.rt = rt
        self.mprop = mprop
        self.rhoprop = rhoprop
        self.FS = FS
        self.Material_NCG = Material_NCG
        self.Material_P = Material_P
        self.Material_V = Material_V
        self.mfrac_tube = mfrac_tube
        self.mfrac_manifold = mfrac_manifold
        self.mfrac_jacket = mfrac_jacket
        self.mfrac_chamber = mfrac_chamber
        self.mfrac_pump = mfrac_pump
        self.mfrac_valve = mfrac_valve
        self.RefMass = RefMass

RL10    = ReferenceEngine(3.20e6, 7.34e4, 61.1,  0.076, 16.85, 351.91, 1.1, Inc_718, D6AC_Steel, D6AC_Steel, 0.0903, 0.2575, 0.0484, 0.2281, 0.3096, 0.0661, 138) #Expander Cycle Reference
LE5     = ReferenceEngine(3.65e6, 1.03e5, 140.0, 0.068, 23.33, 5.5,    1.1, Inc_718, Al_7075_T6, Al_7075_T6, 0.1419, 0.0578, 0.1021, 0.2804, 0.3042, 0.1136, 255)#Gas Generator Cycle Reference #fix rhoprop
SSME    = ReferenceEngine(2.04e7, 2.28e6, 77.5,  0.138, 512.6, 6,      1.2, Inc_718, D6AC_Steel, D6AC_Steel, 0.0907, 0.1984, 0.0737, 0.1955, 0.2763, 0.1654, 3177)#Stage Combustion Cycle Reference #fixrhoprop


#Mass estimation function Nozzle Tubes:
def Mass_Regenerative(Pc, material_N, material_C, material_V, arear, rt, mprop, FS, rhoprop, Cycle): #include Ns and material_P into input if design method is chosen
    #y=

    reference = RL10

    if Cycle == 'EX':
        reference = RL10
    elif Cycle == 'CB':
        reference = RL10
    elif Cycle == 'SC':
        reference = LE5
    else:
        reference = SSME

    #Nozzle Mass
    TubeMass = reference.mfrac_tube*(((Pc/reference.pc)**1)*((material_N.density/reference.Material_NCG.density)**1)*(((material_N.yieldstress_l/FS)/(reference.Material_NCG.yieldstress_l/reference.FS))**(-1))*((arear/reference.arear)**1)*((rt/reference.rt)**2)) #Dimensionless Nozzle Tube Mass
    ManifoldMass = reference.mfrac_manifold*(((Pc/reference.pc)**1)*((material_N.density/reference.Material_NCG.density)**1)*((mprop/reference.mprop)**1)*((material_N.yieldstress_l/reference.Material_NCG.yieldstress_l)**(-1))*((rhoprop/reference.rhoprop)**(-1))*((rt/reference.rt)**1)) #Dimensionless Nozzle Manifold Mass
    JacketMass = reference.mfrac_jacket*(((Pc/reference.pc)**1)*((material_N.density/reference.Material_NCG.density)**1)*(((material_N.yieldstress_l/FS)/(reference.Material_NCG.yieldstress_l/reference.FS))**(-1))*((arear/reference.arear)**1.5)*((rt/reference.rt)**3))#Dimensionless Nozzle Jacket Mass
    
    #Chamber Mass
    ChamberMass = reference.mfrac_chamber*((Pc/reference.pc)**1)*((material_C.density/reference.Material_NCG.density)**1)*(((material_C.yieldstress_l/FS)/(reference.Material_NCG.yieldstress_l/reference.FS))**(-1))*((rt/reference.rt)**1) #Dimensionless Mass of Combustion Chamber + Gas Generator(if applicable)
    
    #TurboPump Mass
    #PumpMass = reference.mfrac_pump*(((Pc/reference.pc)**0.15)*((material_P.density/reference.material_P.density)**1)*(((material_P.yieldstress_l/FS)/(reference.Material_p.yieldstress_l/reference.FS)*(-1)))*((mprop/reference.mprop)**0.9)*((rhoprop/reference.rhoprop)**(-0.45))*((Ns/y)**(-0.6)))
    PumpMass = reference.mfrac_pump*(((Pc/reference.pc)**0.53)*((mprop/reference.mprop)**0.53)*(rhoprop/reference.rhoprop)**(-0.53)) #Dimensionless mass of turbopump (historic data method)

    #Valves
    ValveMass = reference.mfrac_valve*(((Pc/reference.pc)**1.0)*((material_V.density/reference.Material_V.density)**1)*(((material_V.yieldstress_l/FS)/(reference.Material_V.yieldstress_l/reference.FS))**(-1))*((mprop/reference.mprop)**1)*((rhoprop/reference.rhoprop)**(-1)))#Dimensionless mass of Valves
    

    #Total:
    Total_Mass = (TubeMass + ManifoldMass + JacketMass + ChamberMass + PumpMass + ValveMass)*reference.RefMass 

    return Total_Mass


##Cost function 
def Cost(m_engine, R, n):
    TotalCost = 0
    f1  = 1.0
    f2  = 0.27+1.075*R**46.994
    f3  = 1
    f4  = -0.0553*mth.log(n) + 1.0011
    a_m = 4.0

    C_D = 1.1*f1*f2*f3*m_engine**0.58

    F_E = a_m*f4*m_engine**0.46

    TotalCost_MY = C_D + F_E
    TotalCost = TotalCost_MY*200e3
    return TotalCost

# Mass = Mass_Regenerative(3.20e6,Inc_718,Inc_718,D6AC_Steel,61.1,0.076,16.85,1.1,351.91,'EX')
# print(Mass)
