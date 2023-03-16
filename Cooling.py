import math
from scipy.integrate import quad
import scipy.optimize
import scipy.constants
import numpy as np
import Materials
import inspect

# Helper class
class Propellant:
    # Oxidizer
    Ox_name = "LOX"  # Oxidizer name for rocketCEA
    Ox_composition = "O 2"  # Composition of oxidizer for rocketcea
    o_dens = 1141.0  # Oxidizer density
    ocp = 14307.0  # oxidizer cp
    h_ox = -12.979  # oxidizer enthalpy
    o_lamb = 1.0e-3
    omiu = 1.0e-6

    # Fuel
    Fuel_name = "LH2"  # Fuel name for rocketCEA
    Fuel_composition = "H 2"  # Composition of fuel for rocketcea
    f_dens_l = 71.0  # liquid fuel density
    f_dens_g = 1.0  # gaseous fuel density
    f_gamma = 1.4  # fuel gamma
    fcp = 14307.0  # fuel cp
    h_fuel = -9.012  # fuel enthalpy
    R_f = 4.1573  # fuel gas constant
    f_lamb = 1.0e-3
    fmiu = 1.0e-6

    Frozen_state = 0

    # Propellant
    gama = 1.4
    tq = 0.9  # characteristic chemical time of propellant
    MR = 5  # mixture ratio

    def __init__(self, type):  # Placeholder
        if type == 1:
            ox_dens = 1141.0


# To do: calculate the mass flow required in order to end temperature to be equal to operating temperature!


# Holds all the functions to calculate the Heat tranfer in the rocket engine.

# Goal:
#   Calculate the Temperature at the wall, to check if it melts
#   Find the end coolant temperature
#   Find the coolant pressure loss

# Organization:
#   As classes. Each class holds the specific type of cooling (heatsink, radiation,etc...)
#   Each class holds functions to calculate T at the wall:
# Heatsink: Tcalculation()
# RadiationCool: Tcalculation()
# RegenerativeCool: Run()
#   Other functions are defined in each class, which can be used to calculate the heat transmited
#   Each class holds the T at the wall and the heat transmited as a class variable

# All units are in SI


# ------------------Classes-------------------

# General cooling class
# Can be used to organise/hold all types of cooling for a specific application
class Cooling:
    def __init__(self):
        self.heatsink = Heatsink()
        self.radiationcool = RadiationCool()
        self.regencool = RegenerativeCool()


# Heat sink model
# Calculates the Temperature at the wall after a time interval of operation time
# Q is the total energy transfered, and is in J
class Heatsink:
    def __init__(self, Q=-1, m=-1):
        self.Q = Q
        # heat
        self.m = m
        # nozzle mass required
        self.T_calculated = -1

    # Final T after a certain operation time dt, assuming a nozzle mass
    def Tcalculation(self, T0, Tf, h, A, c, dt, m):
        check_positive_args(T0, Tf, h, A, c, dt, m)

        self.T_calculated = (T0 - Tf) * math.e ** (-h * A / (m * c) * dt) + Tf
        check_positive_args(self.T_calculated)

    # Heat absorved assuming a certain nozzle mass
    def Qcalculation(self, T0, Tf, h, A, c, dt, m):
        check_positive_args(T0, Tf, h, A, c, dt, m)

        Q_int_arg = (
            lambda time: h
            * (Tf - ((T0 - Tf) * math.e ** (-h * A / (m * c) * time) + Tf))
            * A
        )
        self.Q = quad(Q_int_arg, 0, dt)
        check_positive_args(self.Q)
    # mass such that the nozzle doesn't melt, assiming Tmelt as the wall temperature
    def mcalculation(self, T0, Tf, h, A, Tmelt, c, dt):
        if Tmelt == Tf:
            raise ValueError(
                "T_melt==Tf, error in mass calculation, for this equality only works after infinite time has passed"
            )
        self.m = h * A * dt / (-math.log((Tmelt - Tf) / (T0 - Tf)) * c)
        check_positive_args(self.m )

# Ratiation cool
# Calculates the equilibrium temperature, taking into account radiated heat
# Q is a power, and is in W
class RadiationCool:
    def __init__(self, Q=-1, t=-1):
        self.Q = Q
        # heat
        self.t = t
        # thickness
        self.T_calculated = -1

    # Calculate the necessary thickness, assuming that the end temperature inside is Tmelt
    def thickcalculation(self, Tmelt, Tr, eps, k, h):
        check_positive_args(Tr, eps, k, Tmelt, h)

        self.t = (
            (h * (Tr - Tmelt) / (eps * scipy.constants.sigma) ** (1 / 4) - Tmelt)
            / (Tr - Tmelt)
            * k
            / h
        )
        self.Q = h(Tr - Tmelt)
        check_positive_args(self.Q,self.t )

    # Defines system of equations required in order to find the end temperature. Auxiliary function
    def Tcalculation_system(self, x, Tr, eps, k, t, h):
        Ti, Tout = x
        return [
            h * (Tr - Ti) - (Tout - Ti) * k / t,
            eps * scipy.constants.sigma * Tout**4 - h * (Tr - Ti),
        ]

    # Calculates the end equilibrium temperature
    # Requires temperature at the wall guesses for the iterative method
    def Tcalculation(self, Toutguess, Tinguess, Tr, eps, k, t):
        check_positive_args(Tr, eps, k, t)

        x0 = [Toutguess, Tinguess]
        sol = scipy.optimize.fsolve(self.Tcalculation_system, x0, args=(Tr, eps, k, t))
        self.T_calculated = sol[0]

        check_positive_args(self.T_calculated)
        return self.T_calculated


# Regenerative Cooling
#
# Main Goal:Get equilibrium wall temperature, end coolant temperature and pressure drop
#
# Options:
# Run(): gives the end temperature of the coolant, the wall temperature and the pressure loss, for a 0D case
# Run1D(): gives the end temperature of the coolant, the wall temperature and the pressure loss, for a 1D case
# Run_for_Toperating0D(): sizes the hydralic diamter to get the wall cooled to operating temperature. Returns the coolant temperature and the pressure loss
# Run_for_Toperating1D(): same as #Run_for_Toperating0D(), but in 1D
# Run1D_iterative_for_m(): Find the massflow required to get the wall cooled to operating temperature. Returns the massflow, the coolant temperature, the wall temperature, and the pressure loss
#
# Note on units
# Q is a power, and is in W
class RegenerativeCool:
    def __init__(self):
        self.Q = 0
        # heat
        self.t = 0
        # thickness
        self.T_calculated = -1

    # ----------------------------------------------------------------------------------------------------------------

    # Calculates and returns equilibrium wall temperature, and end coolant temperature
    # 0D
    # auxiliary function of Run()
    def Tcalculation(self, Tr: float, Ti_co: float, A: float, hg: float):

        q = (Tr - Ti_co) / (1 / hg + self.t / self.Mater.k + 1 / self.hco)
        self.Q += q * A
        Tinext_co = Ti_co + q * A / (self.Prop.fcp * self.m_flow_fuel)
        T_wall = self.t / self.Mater.k * q + Ti_co + q / self.hco

        # Tinext_co: end coolant temperature
        # T_wall: wall temperature
        check_positive_args(Tinext_co, T_wall)
        return Tinext_co, T_wall

    def Tcalculation1D(
        self, Tr: float, Ti_co: float, A: float, hg: float, ArrayCounter: int
    ):
        q = (Tr - Ti_co) / (1 / hg + self.t[ArrayCounter] / self.Mater.k + 1 / self.hco)
        self.Q += q * A[ArrayCounter]
        # print(A[ArrayCounter])
        Tinext_co = Ti_co + q * A[ArrayCounter] / (self.Prop.fcp * self.m_flow_fuel)
        # print("self.hco", self.hco)
        # T_wall = self.t[ArrayCounter] / self.Mater.k * q + Ti_co + q / self.hco
        T_wall = Tr - q / hg
        # Tinext_co: end coolant temperature
        # T_wall: wall temperature
        check_positive_args(Tinext_co, T_wall)
        return Tinext_co, T_wall

    # Calculates pressure loss
    # auxiliary function of Run()
    def pressureloss(self, m_flow_fuel: float, Dr: float, L: float):
        delta_p = self.f * m_flow_fuel**2 / (2 * self.Prop.f_dens_l) * L / Dr
        return delta_p

    def Inicialise(
        self,
        t: float,
        Prop: Propellant,
        Mater: Materials,
        Dr: float,
        Re: float,
        m_flow_fuel: float,
    ):
        self.Q = 0
        # self.Pr = 4 * Prop.f_gamma / (9 * Prop.f_gamma - 5)
        self.Pr = 0.69  # PLACEHOLDER
        self.f = (1.82 * math.log10(Re) - 1.64) ** (-2)
        self.Nu = (
            self.f
            / 8
            * (Re - 1000)
            * self.Pr
            / (1 + 12.7 * math.sqrt(self.f / 8) * (self.Pr**2 / 3 - 1))
        )
        self.hco = self.Nu * Mater.k / Dr

        self.Prop = Prop
        self.m_flow_fuel = m_flow_fuel

        self.t = t
        self.Mater = Mater

    # One of the Main functions for regenerative cooling
    # Takes material properties, flow properties, and coolant pipe properties
    # Returns end wall temperature, end coolant temperature, pressure drop
    # Saves/Updates Q
    def Run(
        self,
        Tr: np.array,
        hg: float,
        t: float,
        Prop: Propellant,
        Mater: Materials,
        Dr: float,
        A: float,
        Ti_co: float,
        Re: float,
        m_flow_fuel: float,
        L: float,
    ):
        check_positive_args(Tr, hg, t, Dr, A, Ti_co, m_flow_fuel, Re, L)

        self.Inicialise(t, Prop, Mater, Dr, Re, m_flow_fuel)

        T_co_calcualted, T_wall_calcualted = self.Tcalculation(Tr, Ti_co, A, hg)
        ploss = self.pressureloss(m_flow_fuel, Dr, L)

        # T_co_list=[0 for i in range(len(Tr))]
        # T_co_list[0]=Ti_co
        # for i in Tr:
        #   T_co_calcualted[i], T_wall_calcualted[i] = self.Tcalculation(Tr[i],T_co_list[i],A[i])
        check_positive_args(T_co_calcualted, T_wall_calcualted, ploss)
        return T_co_calcualted, T_wall_calcualted, ploss

    # Run function, but for 1D case
    def Run1D(
        self,
        Tr: np.array,
        hg: float,
        t: float,
        Prop: Propellant,
        Mater: Materials,
        Dr,
        A,
        Ti_co,
        Re,
        m_flow_fuel,
        L,
    ):
        check_positive_args(Tr, hg, t, Dr, A, Ti_co, m_flow_fuel, Re, L)

        # Inicialise variables
        self.Inicialise(t, Prop, Mater, Dr, Re, m_flow_fuel)
        self.D = Dr

        # T_co_calcualted, T_wall_calcualted = self.Tcalculation(Tr, Ti_co, A, hg)

        # Calculate pressure loss
        ploss = self.pressureloss(m_flow_fuel, Dr, L)

        # Inicialise variables with placeholders
        T_co_calcualted = [0 for i in range(len(Tr) + 1)]
        T_co_calcualted[0] = Ti_co
        T_wall_calcualted = [0 for i in range(len(Tr))]

        # Calculate the wall temperature and the coolant temperature for each point along the wall
        for i in range(len(Tr)):
            T_co_calcualted[i + 1], T_wall_calcualted[i] = self.Tcalculation1D(
                Tr[i], T_co_calcualted[i], A, hg[i], i
            )
        # print(T_co_calcualted)
        check_positive_args(T_co_calcualted, T_wall_calcualted, ploss)
        return T_co_calcualted, T_wall_calcualted, ploss

    # ----------------------------------------------------------------------------------------------------------------

    # One of the Main functions for regenerative cooling
    # These solvers size the hydralic diameter such that it runs at operating temperature
    # This is for only one tube along one wall, mass flow might have to be adjusted
    def Run_for_Toperating0D(
        self,
        Tr: np.array,
        hg: float,
        t: float,
        Prop: Propellant,
        Mater: Materials,
        A,
        Ti_co,
        m_flow_fuel,
        L,
    ):
        check_positive_args(Tr, hg, t, A, Ti_co, m_flow_fuel, L)

        # Inicicialise variables
        self.Prop = Prop
        self.m_flow_fuel = m_flow_fuel
        self.t = t
        self.Mater = Mater

        Twh = self.Mater.OpTemp_u
        if Tr < Twh:
            print(
                "Temperature at the wall is smaller than operating temperature. No need for regenerative cooling"
            )
            self.D = -404
            return Ti_co, 0
        if Ti_co > Tr:
            raise Exception("Ti_co > Tr, Ti_co: ", Ti_co)

        # calculate the convection coefficient for the coolant
        self.hco = (Tr - Twh) / (
            (Twh - Ti_co) * (1 / hg + t / self.Mater.k)
            - (Tr - Ti_co) * t / self.Mater.k
        )

        # print("data: ", self.hco)

        # Find the hydralic diameter which gives the correct convection coefficient for the coolant
        D0 = 0.00001
        D = scipy.optimize.fsolve(self.SolveForD, D0)

        # Update coolant temperature
        q = (Tr - Ti_co) / (1 / hg + self.t / self.Mater.k + 1 / self.hco)
        self.Q += q * A
        T_co_calcualted = Ti_co + q * A / (self.Prop.fcp * self.m_flow_fuel)
        # print("h_co: ", self.hco)
        # T_co_calcualted = Ti_co + q * A / (2*10**6* self.m_flow_fuel)

        # Calculate Pressure loss
        ploss = self.pressureloss(m_flow_fuel, D, L)

        self.D = D
        check_positive_args(T_co_calcualted, ploss)
        return T_co_calcualted, ploss

    # Find the hydralic diameter which gives the correct convection coefficient for the coolant
    def SolveForD(self, D: float):

        # Calculate the required parameters to calculate the h_co
        Re = self.m_flow_fuel / self.Prop.fmiu * 4 / (math.pi * D)
        self.Pr = 1  # PLACEHOLDER
        self.f = (1.82 * math.log10(Re) - 1.64) ** (-2)
        self.Nu = (
            self.f
            / 8
            * (Re - 1000)
            * self.Pr
            / (1 + 12.7 * math.sqrt(self.f / 8) * (self.Pr ** (2 / 3) - 1))
        )

        # Compare the calculate h_co with the correct one
        eq = self.Nu * self.Mater.k / D - self.hco
        return eq

    # One of the Main functions for regenerative cooling
    def Run_for_Toperating1D(
        self,
        Tr_array: np.array,
        hg: float,
        t: float,
        Prop: Propellant,
        Mater: Materials,
        A,
        Ti_co,
        m_flow_fuel,
        L,
        y,
    ):
        check_positive_args(Tr_array, hg, t, A, Ti_co, m_flow_fuel, L)

        self.Prop = Prop
        self.m_flow_fuel = m_flow_fuel
        self.t = t
        self.Mater = Mater

        if self.Mater.OpTemp_u > float(np.amax(Tr_array)):
            print(
                "Temperature at the wall is smaller than operating temperature. No need for regenerative cooling"
            )
            self.D = -404
            return Ti_co, 0

        zeroDcool = RegenerativeCool()
        Ti_co_array = [0 for i in range(len(Tr_array) + 1)]
        ploss = [0 for i in range(len(Tr_array))]
        Ti_co_array[0] = Ti_co
        D = 100000000000000000000

        # Iterative loop
        # Find the coolant temperature and the pressure loss along the wall
        # optimise for D, to get operating temperature
        # if the a new calculated D is smaller than the old one, update all values from the begining before continuing
        for i in range(len(Tr_array)):
            # update the surface area for each node along the wall, for each iteration
            A = (2 * math.pi * y[i]) * L / len(Tr_array)

            if D > 0.1:
                D = 0.1

            # optimise for D, to get operating temperature
            Ti_co_array[i + 1], ploss[i] = zeroDcool.Run_for_Toperating0D(
                Tr_array[i],
                hg[i],
                t[i],
                Prop,
                Mater,
                A,
                Ti_co_array[i],
                m_flow_fuel,
                L / len(Tr_array),
            )
            if zeroDcool.D < D:
                D = zeroDcool.D
                # A = (2 * math.pi * y[i]) * L / len(Tr_array)
                zeroDcool.Q = 0
                if D > 0.1:
                    D = 0.1
                for j in range(i + 1):
                    Ti_co_array[j + 1], ploss[j] = zeroDcool.Run_for_Toperating0D(
                        Tr_array[j],
                        hg[j],
                        t[j],
                        Prop,
                        Mater,
                        A,
                        Ti_co_array[j],
                        m_flow_fuel,
                        L / len(Tr_array),
                    )
        self.D = D
        self.Q = zeroDcool.Q
        # print(Ti_co_array)
        check_positive_args(Ti_co_array[len(Ti_co_array) - 1], ploss)
        return Ti_co_array[len(Ti_co_array) - 1], ploss

    # ----------------------------------------------------------------------------------------------------------------

    def R1D_findm_for_TOp(
        self,
        m_flow_fuel: float,
        Tr: float,
        hg: float,
        t: float,
        Prop: Propellant,
        Mater: Materials,
        Dr: float,
        A: float,
        Ti_co: float,
        Re: float,
        L: float,
    ):
        check_positive_args(Tr, hg, t, Dr, A, Ti_co, Re, L)

        # Inicialise variables, namely those required for heat transfer, such as Nu, Pr, hco,...
        self.Inicialise(t, Prop, Mater, Dr, Re, m_flow_fuel)
        self.D = Dr
        # T_co_calcualted, T_wall_calcualted = self.Tcalculation(Tr, Ti_co, A, hg)
        ploss = self.pressureloss(m_flow_fuel, Dr, L)

        T_co_calcualted = [0 for i in range(len(Tr) + 1)]
        T_co_calcualted[0] = Ti_co
        T_wall_calculated = [0 for i in range(len(Tr))]

        # Calculate the wall temperature and the coolant temperature for each point along the wall
        for i in range(len(Tr)):
            T_co_calcualted[i + 1], T_wall_calculated[i] = self.Tcalculation1D(
                Tr[i], T_co_calcualted[i], A, hg[i], i
            )
        # print(T_co_calcualted)
        self.T_col = T_co_calcualted
        self.Tw_wall_calculated = T_wall_calculated
        return max(T_wall_calculated) - Mater.OpTemp_u

    # One of the Main functions for regenerative cooling
    # Finds the mass flow for the wall to be at operation temperature
    # returns the mass flow, the coolant temperature, the wall temperature, and the pressure loss ( ploss as a float)
    def Run1D_iterative_for_m(
        self,
        Tr: np.array,
        hg: float,
        t: float,
        Prop: Propellant,
        Mater: Materials,
        Dr: float,
        A: float,
        Ti_co: float,
        Re: float,
        L: float,
    ):
        check_positive_args(Tr, hg, t, Dr, A, Ti_co, Re, L)

        # Optimise mass flow for operating wall temperature
        m_flow_fuel = scipy.optimize.fsolve(
            self.R1D_findm_for_TOp,
            1,
            args=(Tr, hg, t, Prop, Mater, Dr, A, Ti_co, Re, L),
        )

        # Calculate pressure loss; ploss is a float
        ploss = self.pressureloss(m_flow_fuel, Dr, L)

        check_positive_args(m_flow_fuel, self.T_col, self.Tw_wall_calculated, ploss)

        return m_flow_fuel, self.T_col, self.Tw_wall_calculated, ploss


# Sanitisaion function


def check_positive_args(*args):
    for arg in args:
        if isinstance(arg, (int, float, np.int32, np.generic)):
            if arg < 0:
                raise ValueError("All numerical arguments must be positive")
        elif isinstance(arg, (list, tuple, np.ndarray)):
            if any(x < 0 for x in arg):
                raise ValueError("All elements of numerical arguments must be positive")
        else:
            raise ValueError("Unsupported argument type: {type(arg).__name__}")
