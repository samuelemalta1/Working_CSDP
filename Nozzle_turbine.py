def Turbine_nozzle(m_p,Pc,Prop,Pamb,Default,h_fuel,h_ox,rho_fuel,rho_ox,MR):

    Ox=Prop.Ox_name
    Fuel=Prop.Fuel_name
    frozen_state=Prop.Frozen_state
    De_turbine_noz_max=Default.De_turbine_noz_max
    Ox_composition=Prop.Ox_composition
    Fuel_composition=Prop.Fuel_composition
    Pc = Pc/1.0e5
    Pamb = Pamb/1.0e5
    

    from rocketcea.cea_obj_w_units import CEA_Obj
    from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer
    import math as mth
    import matplotlib.pyplot as plt
    import numpy as geek

    card_str="""
    oxid {} {} 
    h,kj/kg={} t(k)={} rho={}
    """.format(Ox,Ox_composition,h_ox,298,rho_ox)
    add_new_oxidizer('Ox_turbine',card_str)

    card_str="""
    fuel {} {} 
    h,kj/kg={} t(k)={} rho={}
    """.format(Fuel,Fuel_composition,h_fuel,298,rho_fuel)

    add_new_fuel('Fuel_turbine',card_str)
    
    ispObj = CEA_Obj( oxName='Ox_turbine', fuelName='Fuel_turbine',cstar_units='m/s',pressure_units='bar',temperature_units='K',isp_units='sec',density_units='kg/m^3',specific_heat_units='J/kg-K',viscosity_units='poise')

    Ae_max=De_turbine_noz_max**2.0/4.0*mth.pi
    At=ispObj.get_Cstar(Pc,MR)*m_p/(Pc*1.0e5)
    variation=0.1
    v_eff=-1.0
    it=0

    while v_eff<0.0 or variation>0.01:

        eps_max=Ae_max/At # Maximum expansion ratio with this throat area
        Pratio_max=ispObj.get_PcOvPe(Pc=Pc,MR=MR,eps=eps_max,frozen=frozen_state,frozenAtThroat=frozen_state)
        Pe_max=Pc/Pratio_max

        if Pe_max>Pamb:
            difference=0.001
            Pe=Pe_max
            Ae=Ae_max; #In this case the maximum exit area cannot reach adapted conditions and we will simply take that
        else:
            difference=0.1
            Ae_1=At+0.000001; # In this case we enter the loop and start the iteration with Minimum exit area with this throat area
        
        while difference>0.01:

            eps_1= Ae_1/At # First expansion ratio 
            Pratio_1=ispObj.get_PcOvPe(Pc=Pc,MR=MR,eps=eps_1,frozen=frozen_state,frozenAtThroat=frozen_state)
            Pe_1=Pc/Pratio_1 # Exit pressure 1
            
            difference=abs(Pe_1-Pamb)/Pamb 

            if difference>0.01:

                Ae_2=Ae_1*(Pe_1/Pamb)**0.3 # Change value for next iteration
                Ae_1=Ae_2;
            else:
                Ae=Ae_1 # Finish iteration if we have reached the convergence
                Pe=Pe_1;
        eps_actual=Ae/At;
        
        Isp_it=ispObj.estimate_Ambient_Isp(Pc=Pc,MR=MR,eps=eps_actual,Pamb=Pamb,frozen=frozen_state,frozenAtThroat=frozen_state) # Calculates Isp for this iteration
        v_eff_it=Isp_it[0]*9.80665; # Calculates effective velocity for this iteration

        sonic_vel_it=ispObj.get_SonicVelocities(Pc=Pc,MR=MR,eps=eps_actual,frozen=frozen_state,frozenAtThroat=frozen_state)
        u_t=sonic_vel_it[1] # Velocity of sound in the throat

        rhos=ispObj.get_Densities(Pc=Pc,MR=MR,eps=eps_actual,frozen=frozen_state,frozenAtThroat=frozen_state)
        rho_t=rhos[1]  # Density in the throat

        m_p_it=rho_t*u_t*At # Mass flow rate from continuity equation in the throat
        variation=abs(m_p_it-m_p)/m_p # Variation from target force

        if variation>0.01: # If we have to iterate again, we change the throat area
            At_it=At
            At=At_it*m_p/m_p_it;
    
        it=it+1;
        if it>10000:
            break;
        v_eff=v_eff_it;
    
    Isp=Isp_it[0]
    return Isp
