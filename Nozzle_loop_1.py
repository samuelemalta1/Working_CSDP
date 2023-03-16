def Nozzle_loop_1(Pc,F_tar,Pamb,Propellant,Default,Nozzle_type):
    ## INPUTS:
    # Pc= Chamber pressure in bar
    # F=Thrust in Newtons
    # Pamb= Ambient pressure in bar
    # Propellant.Ox_name= string with oxydizer name for RocketCEA
    # Propellant.Fuel_name= string with fuel name for Rocket CEA
    # Propellant.Frozen_state: 0= equilibrium everywhere, 1= equilibrium in the convergent and frozen in the divergent
    # Default.MR: 0=calculate for best c*, 'every_other_value'=Mixture ratio to use
    # Default.De_max= maximum exit diameter in meters
    # Default.toll_c_star= tollerance for the best MR calculation
    # Default.toll_F_obj= tollerance on the target force to be achieved during mass flow rate calculation
    # Default.Max_iterations_mass_flow= max iterations for the throat area and mass flow calculations
    # Default.toll_P_adapted= tollerance on difference between exit pressure and ambient pressure if nozzle can be adapted
    # Nozzle_type= type of nozzle, 0== Conical, 1== Bell

    ## OUTPUTS:
    # m_p = mass flow rate kg/s
    # Tc= combustion chamber temperature K
    # MR= mixture ratio
    # At= throat area m^2
    # eps= expansion ratio
    # Isp= Specific impulse s

    from rocketcea.cea_obj_w_units import CEA_Obj
    import math as mth
    import matplotlib.pyplot as plt
    import numpy as geek
    
    Ox=Propellant.Ox_name
    Fuel=Propellant.Fuel_name
    frozen_state=Propellant.Frozen_state

    MR=Default.MR
    De_max=Default.De_max
    toll_c_star=Default.toll_c_star
    toll_F_obj=Default.toll_F_obj
    Max_iterations_mass_flow=Default.Max_iterations_mass_flow
    toll_P_adapted=Default.toll_P_adapted
    theta_conical=mth.radians(Default.Theta_conical)
    Theta_bell=mth.radians(Default.Theta_bell)
    TH_exit_bell=mth.radians(Default.TH_exit_bell)
    Ru_bell=Default.R_u_bell
    #eps_m=Default.eps_m
    eps_m=77.5

    ispObj = CEA_Obj( oxName=Ox, fuelName=Fuel,cstar_units='m/s',pressure_units='bar',temperature_units='K',isp_units='sec',density_units='kg/m^3',specific_heat_units='J/kg-K',viscosity_units='poise',thermal_cond_units='W/cm-degC')


    # Sanitizing the inputs

    error=0;

    # Target Thrust
    if F_tar>1000000:
        print("Warning, the target thrust is higher than the range of thrusts this tool is intended for")
    #if F<=0:
    #    print("ERROR: Negative thrust")
    #error=1;
    #    return 0,0,0,0,0,0,0,0,0,0,0,error
    
    # Chamber pressure
    if Pc>300:
        print("Warning, the chamber pressure is very high, it might not be feasible")
    if Pc<=Pamb:
        print("ERROR, Chamber pressure lower than ambient pressure")
        error=1
        return 0,0,0,0,0,0,0,0,0,0,0,error
    
    
    # Ambient pressure
    #if Pamb>1.01325:
        print("Warning, ambient pressure higher than sea level standard pressure");
    
    # MR
    #if MR<0:
    #    print("ERROR, mass flow rate is lower than zero")
    #    error=1
    #    return 0,0,0,0,0,0,0,0,0,0,0,error
    
    #if MR>20:
    #    print("Warning, mass flow rate is very high, program might crash due to unphysical conditions")
    

    # Losses due to divergent part in conical nozzle
    if Nozzle_type==0:
        eps_loss=0.5*(1-mth.cos(theta_conical))
        F=F_tar/(1-eps_loss);
    
    if MR==0:
        MR_1=0.01
        MR_2=150

        diff=MR_1-MR_2
        MR_prev=MR_1;

        while abs(diff)>toll_c_star:
            
            MR_curr=(MR_2+MR_1)/2
            diff=MR_curr-MR_prev

            c_star_curr=ispObj.get_Cstar(Pc=Pc,MR=MR_curr)

            MR_curr_low=MR_curr-0.001
            MR_curr_high=MR_curr+0.001

            c_star_curr_low=ispObj.get_Cstar(Pc=Pc,MR=MR_curr_low)

            if c_star_curr>c_star_curr_low:
                MR_1=MR_curr;
            else:
                MR_2=MR_curr;
            MR_prev=MR_curr;
        MR=MR_curr #Mixture ratio
        c_star=c_star_curr;
    else:
        c_star=ispObj.get_Cstar(Pc=Pc,MR=MR);

    Tc=ispObj.get_Tcomb(Pc=Pc,MR=MR) # Function that returns the combustion chamber temperature

    v_eff=-1
    it=0

    At=0.001 # Initial value for throat area, to be better defined !!!!!!!!!!!!!!

   
    Ae_max=mth.pi*De_max**2/4 #Maximum exit area of the nozzle

 
    variation=10 # Value to enter the loop

    while v_eff<0 or variation>toll_F_obj:

        eps_max=Ae_max/At # Maximum expansion ratio with this throat area
        if eps_max>eps_m:
            eps_max=eps_m;
        
        Pratio_max=ispObj.get_PcOvPe(Pc=Pc,MR=MR,eps=eps_max,frozen=frozen_state,frozenAtThroat=frozen_state)
        Pe_max=Pc/Pratio_max

        if Pe_max>Pamb or eps_max==eps_m:
            difference=0.001
            Pe=Pe_max
            Ae=At*eps_max; #In this case the maximum exit area cannot reach adapted conditions and we will simply take that
        else:
            difference=0.1
            Ae_1=At+0.001; # In this case we enter the loop and start the iteration with Minimum exit area with this throat area
        
        while difference>toll_P_adapted:

            eps_1= Ae_1/At # First expansion ratio 
            
            Pratio_1=ispObj.get_PcOvPe(Pc=Pc,MR=MR,eps=eps_1,frozen=frozen_state,frozenAtThroat=frozen_state)
            Pe_1=Pc/Pratio_1 # Exit pressure 1
            
            difference=abs(Pe_1-Pamb)/Pamb 

            if difference>0.01:
                
                #g_it=ispObj.get_Chamber_MolWt_gamma(Pc=Pc,MR=MR,eps=eps_1)
                #g_cur=g_it[1]
                #G_cur=mth.sqrt(g_cur)*((2/(g_cur+1))**((g_cur+1)/2*(g_cur-1)))
                #eps_it=G_cur/(mth.sqrt(2*g_cur/(g_cur-1)*(Pamb/Pc)**(2/g_cur)*(1-(Pamb/Pc)**((g-1)/g))))
                #Ae_1=At*eps_it;
                Ae_2=Ae_1*Pe_1/Pamb # Change value for next iteration
                Ae_1=Ae_2;
            else:
                Ae=Ae_1 # Finish iteration if we have reached the convergence
                Pe=Pe_1;
            
        eps_actual=Ae/At;
        
        # Divergence losses in case of a bell nozzle:
        if Nozzle_type==1:
            R_t=mth.sqrt(At/mth.pi)
            ye=mth.sqrt(eps_actual)*R_t
            xp=Ru_bell*R_t*mth.sin(Theta_bell)
            yp=(1+Ru_bell)*R_t-Ru_bell*R_t*mth.cos(Theta_bell)
            a=(mth.tan(mth.pi/2-TH_exit_bell)-mth.tan(mth.pi/2-Theta_bell))/(2*(ye-yp))
            b=mth.tan(mth.pi/2-Theta_bell)-2*(mth.tan(mth.pi/2-TH_exit_bell)-mth.tan(mth.pi/2-Theta_bell))/(2*(ye-yp))*yp
            c=xp-a*yp**2-b*yp
            L_nozzle_div=a*ye**2+b*ye+c
            alpha=mth.atan((ye-yp)/L_nozzle_div)
            eps_loss=0.5*(1-mth.cos((alpha+TH_exit_bell)/2))
            F=F_tar/(1-eps_loss);
        
        Isp_it=ispObj.estimate_Ambient_Isp(Pc=Pc,MR=MR,eps=eps_actual,Pamb=Pamb,frozen=frozen_state,frozenAtThroat=frozen_state) # Calculates Isp for this iteration
        v_eff_it=Isp_it[0]*9.80665*(1-eps_loss); # Calculates effective velocity for this iteration

        m_p_it=Pc*100000*At/c_star # Mass flow rate from continuity equation in the throat

        F_it=m_p_it*v_eff_it # Force computed at this iteration
    
        variation=abs(F_it-F)/F # Variation from target force

        if variation>toll_F_obj: # If we have to iterate again, we change the throat area
            mp2=m_p_it*(F/F_it);
            At=c_star*mp2/(Pc*100000);
    
        it=it+1;
        if it>Max_iterations_mass_flow:
            break;
        v_eff=v_eff_it;
    
    m_p=m_p_it # Propellant mass flow
    eps=Ae/At # Expansion ratio
    Isp=Isp_it # Isp value

    # Combustion chamber properties

    rhos=ispObj.get_Densities(Pc=Pc,MR=MR,eps=eps_actual,frozen=frozen_state,frozenAtThroat=frozen_state)
    Transp_c=ispObj.get_Chamber_Transport(Pc=Pc,MR=MR,eps=eps_actual,frozen=frozen_state)

    rho_c=rhos[0]
    cp_c=Transp_c[0]
    mu_c=Transp_c[1]
    k_c=Transp_c[2]/100
    Pr_c=Transp_c[3]

    ## Sanitizing outputs

    # Mass flow rate
    if m_p<=0:
        print("ERROR: mass flow rate is lower than 0")
        error=1

    
    # Throat area
    if At<=0:
        print("ERROR, throat area lower than 0")
        error=1

    
    # Expansion ratio
    if eps<=1:
        print("Error, expansion ratio is lower than 1")
        error=1
        
    
    # Isp
    if Isp[0]<=0:
        print("ERROR, Isp<0")
        error=1
        
    
    # chamber properties

    if rho_c<=0 or cp_c<=0 or mu_c <=0 or k_c<=0 or Pr_c<=0:
        print("ERROR, chamber properties have negative value")
        error=1
        
    
    return m_p,Tc,MR,At,eps,Isp[0]*(1-eps_loss),rho_c,cp_c,mu_c,k_c,Pr_c,error



