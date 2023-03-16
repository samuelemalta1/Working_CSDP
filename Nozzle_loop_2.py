def Nozzle_loop(Pc,Tc,Propellant,Material,Nozzle_type,MR,eps,At,m_p,Dc,Default):
    ## INPUTS:
    # Default
    # Pc
    # Tc
    # Propellant
    # Material
    # Nozzle_type
    # MR
    # eps
    # At
    # m_p
    # Dc
    ## Outputs
    # t_noz Thickness of the nozzle
    # x_noz Points on the axis at which properties are evaluated
    # y_noz Geometry of the nozzle
    # Tw_ad_noz Adiabatic wall temperature in the nozzle
    # h_c_noz Convective heat coefficient in the nozzle
    # P_noz Pressure in the nozzle
    # T_noz Temperature in the nozzle

    from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer
    from rocketcea.cea_obj_w_units import CEA_Obj
    import math as mth
    import matplotlib.pyplot as plt
    import numpy as geek
    
    # We begin by importing all the defaul values and the propellant properties, as well as setting the CEA environment
    Ox=Propellant.Ox_name
    Fuel=Propellant.Fuel_name
    frozen_state=Propellant.Frozen_state
    
    ispObj = CEA_Obj( oxName=Ox, fuelName=Fuel,useFastLookup=0, makeOutput=0,isp_units='sec',cstar_units='m/sec',pressure_units='Bar',temperature_units='K', sonic_velocity_units='m/s',enthalpy_units='kJ/kg',density_units='kg/m^3',specific_heat_units='J/kg-K',viscosity_units='poise',thermal_cond_units='W/cm-degC',fac_CR=None, make_debug_prints=False)

    ## Sanitizing inputs
    error=0;
    # Angles
    if Default.Theta_con>=90:
        print("ERROR: Convergent angle higher than 90 degrees")
        error=1
        return 0,0,0,0,0,0,0,0,0,0,0,0,error;

    if Default.Theta_con>70:
        print("Warning: Convergent angle very high, might cause real life issues")
    
    if Default.Theta_conical>20:
        print("Warning: Divergent angle too high, might cause flow separation")
    
    if Default.Theta_conical>=90:
        print("ERROR: Divergent angle higher than 90 degrees")
        error=1
        return 0,0,0,0,0,0,0,0,0,0,0,0,error;

    if Default.Theta_bell>65:
        print("Warning: Divergent angle very high, might cause real life issues")
    
    if Default.Theta_bell>=90:
        print("ERROR: Divergent angle higher than 90 degrees")
        error=1
        return 0,0,0,0,0,0,0,0,0,0,0,0,error;

    if Default.TH_exit_bell>15:
        print("Warning: exit angle is very high, might lead to high losses")
    
    if Default.TH_exit_bell>Default.Theta_bell:
        print("ERROR: Exit angle higher than divergent angle")
        error=1
        return 0,0,0,0,0,0,0,0,0,0,0,0,error;

    # Radius of curvature

    if Default.R_u_bell > 1.5:
        print("Warning: throat curvature angle very high")
    
    if Default.R_u_bell<=0:
        print("ERROR: throat curvature angle lower than 0")
        error=1
        return 0,0,0,0,0,0,0,0,0,0,0,0,error;

    if Default.R_u_bell<0.5:
        print("Warning: throat curvature angle very low, might lead to flow separation")

    # Chamber diameter
    if (mth.pi*Dc**2/4)<At:
        print("ERROR: chamber diameter smaller than throat diameter")
        error=1
        return 0,0,0,0,0,0,0,0,0,0,0,0,error;

    # Nozzle resolution
    if Default.noz_res<100:
        print("Warning: nozzle resolution is low, might lead to inaccurate results")
    
    if Default.noz_res<10:
        print("ERROR: nozzle resolution is too low")
        error=1
        return 0,0,0,0,0,0,0,0,0,0,0,0,error;

    if Default.n_cool>99:
        print("Warning, too many points for cooling properties resolution can slow down the code")
    
    if Default.n_cool<10:
        print("ERROR: nozzle cooling resolution to low")
        error=1
        return 0,0,0,0,0,0,0,0,0,0,0,0,error;
    
    T_w=Material.OpTemp_u
    theta_con=mth.radians(Default.Theta_con)
    Theta_conical=mth.radians(Default.Theta_conical)
    Theta_bell=mth.radians(Default.Theta_bell)
    TH_exit_bell=mth.radians(Default.TH_exit_bell)
    Ru_bell=Default.R_u_bell
    num_cool=Default.n_cool

    noz_res=Default.noz_res
    # Definition of the geometry of the nozzle
    R_t=mth.sqrt(At/(mth.pi)) #Throat radius
    R_u=R_t*Default.R_u_ratio 
    con_ratio=Dc**2/((2*R_t)**2) # Contraction ratio (from chamber to throat)

    ## Subdivision of the discretization between the convergent, throat and divergent parts of the nozzle
    num_con=round(0.15*noz_res)
    num_th=round(0.25*noz_res)
    num_div=round(0.6*noz_res)

    Dt=R_t*2; # Throat diameter
    De=Dt*(mth.sqrt(eps)) # Exit diameter
    if Nozzle_type==0: # If to determine whether the nozzle is conical (==0) or bell (==1)
        L_nozzle_div=(mth.sqrt(eps)*R_t-R_t-R_u*(1-mth.cos(mth.cos(Theta_conical))))/(mth.tan(Theta_conical))+R_u*mth.sin(Theta_conical) # Length of the divergent part of the nozzle for the conical nozzle
        L_nozzle_con=(mth.sqrt(con_ratio)*R_t-R_t-R_u*(1-mth.cos(mth.cos(theta_con))))/(mth.tan(theta_con))+R_u*mth.sin(theta_con) # Length of the convergent part of the nozzle for the conical
        L_tot=L_nozzle_con+L_nozzle_div # Total length of the nozzle

        xp=L_nozzle_con+R_u*mth.sin(Theta_conical) # point P as defined in the documentation (x value)
        yp=R_t+(1-mth.cos(Theta_conical))*R_u # point P as defined in the documentation (y value)

        xp_con=L_nozzle_con-R_u*mth.sin(theta_con) # point P as defined in the documentation, but for the convergent(x value)
        yp_con=R_t+(1-mth.cos(theta_con))*R_u # point P as defined in the documentation, but for the convergent(y value)
        
        n=num_con 
        x_con=geek.linspace(0,xp_con,num=n) # Discretization of points in the convergent
        a_con=(yp_con-Dc/2)/(xp_con) # Coefficient for the convergent part
        y_con=Dc/2+a_con*x_con # With this part we have defined the coordinates for the convergent geometry

        # Discretization of the first part of the throat (convergent section of the throat area)
        n1=round(num_th/2)
        x_throat1=geek.linspace(xp_con,L_nozzle_con,num=n1)
        th_step=[]
        y_throat1=[]
        for i in x_throat1:
            th_step_cur=mth.asin((L_nozzle_con-i)/R_u)
            th_step.append(th_step_cur)
            y_cur=R_t+(1-mth.cos(th_step_cur))*R_u
            y_throat1.append(y_cur) # With this part we have defined the coordinates for the first part of the throat using the radius of curvature

        # Same process of discretizaion but for the divergent section of the throat
        n2=round(num_th/2)
        x_throat2=geek.linspace(L_nozzle_con,xp,num=n2)#With this part we have defined the coordinates for the second part of the throat

        th_step=[]
        y_throat2=[]
        for i in x_throat2:
            th_step_cur=mth.asin((i-L_nozzle_con)/R_u)
            th_step.append(th_step_cur)
            y_cur=R_t+(1-mth.cos(th_step_cur))*R_u
            y_throat2.append(y_cur)
        
        # Discretization of the divergent part of the nozzle. Same processes as for the convergent, following documentation.
        n3=num_div
        x_div=geek.linspace(xp,L_tot,num=n3)
        a_div=(mth.sqrt(eps)*R_t-yp)/(L_tot-xp)
        y_div=yp+a_div*(x_div-xp) # With this part we have defined the coordinates for the divergent part of the nozzle
        x_2=geek.concatenate((x_throat2,x_div))

        ## Limitation of the points in the divergent part to avoid reading CEA file too many times
        if (n2+n3)>num_cool: # Check whether the number of points exceeds the limit set
            e_tran_Num=num_cool
            num_th_2=round(len(x_throat2)/len(x_2))*e_tran_Num
            #A_div_cool=geek.linspace(A_div[0],A_div[-1],e_tran_Num)
            x_th2_cool=geek.linspace(x_throat2[0],x_throat2[-1],num_th_2)
            x_div_cool=geek.linspace(x_div[0],x_div[-1],(e_tran_Num-num_th_2)) #Re-discretization of the divergent for cooling properties
            y_th2_cool=[]
            for i in x_th2_cool:
                th_step_cur=mth.asin((i-L_nozzle_con)/R_u)
                th_step.append(th_step_cur)
                y_cur=R_t+(1-mth.cos(th_step_cur))*R_u
                y_th2_cool.append(y_cur);
            a_div_cool=(mth.sqrt(eps)*R_t-yp)/(L_tot-xp)
            y_div_cool=yp+a_div_cool*(x_div_cool-xp);
        else:
            x_div_cool=geek.zeros(len(x_div)) # If the number of points in the divergent does not exceed the limit, we keep the previous discretization
            for i in range(len(x_div_cool)):
                x_div_cool[i]=x_div[i];
            x_th2_cool=geek.zeros(len(x_th2_cool))
            for i in range(len(x_th2_cool)):
                x_th2_cool[i]=x_throat2[i];
            y_th2_cool=geek.zeros(len(y_throat2))
            for i in range(len(y_th2_cool)):
                y_th2_cool[i]=y_throat2[i];
            y_div_cool=geek.zeros(len(y_div))
            for i in range(len(y_div_cool)):
                y_div_cool[i]=y_div[i];
        
        # Coordinates of the discretizaion of the nozzle (no limits for cooling)
        x_noz=geek.concatenate((x_con,x_throat1,x_throat2,x_div))
        y_noz=geek.concatenate((y_con,y_throat1,y_throat2,y_div)) 

        # Coordinate of the discretizaion of the nozzle (limits for the cooling)
        x_noz_cool=geek.concatenate((x_con,x_throat1,x_th2_cool,x_div_cool))
        y_noz_cool=geek.concatenate((y_con,y_throat1,y_th2_cool,y_div_cool))
    else: ## Definition of bell nozzle geometry
        ye=mth.sqrt(eps)*R_t # y coordinate for the exit section of the nozzle
        xp=Ru_bell*R_t*mth.sin(Theta_bell) #Point P as defined in the documentation (x value)
        yp=(1+Ru_bell)*R_t-Ru_bell*R_t*mth.cos(Theta_bell)#Point P as defined in the documentation (y value)
        a=(mth.tan(mth.pi/2-TH_exit_bell)-mth.tan(mth.pi/2-Theta_bell))/(2*(ye-yp)) # Definition of parabolic coefficients
        b=mth.tan(mth.pi/2-Theta_bell)-2*(mth.tan(mth.pi/2-TH_exit_bell)-mth.tan(mth.pi/2-Theta_bell))/(2*(ye-yp))*yp
        c=xp-a*yp**2-b*yp
        L_nozzle_div=a*ye**2+b*ye+c # Length of divergent part of the nozzle
        L_nozzle_con=((mth.sqrt(con_ratio)-1)*R_t+R_u*(1/mth.cos(theta_con)-1))/mth.tan(theta_con) # Length of convergent part of the nozzle
        L_tot=L_nozzle_con+L_nozzle_div #Total length of the nozzle
        
        ## Same discretization process as for the convergent, with same check for discretization points exceeding the cooling limit
        xp=xp+L_nozzle_con

        xp_con=L_nozzle_con-R_u*mth.sin(theta_con)
        yp_con=R_t+(1-mth.cos(theta_con))*R_u
        
        n=num_con
        x_con=geek.linspace(0,xp_con,num=n)
        a_con=(yp_con-Dc/2)/(xp_con)
        y_con=Dc/2+a_con*x_con # With this we have defined the geometry of the convergent part

        y_throat1=[]
        th_step=[]

        n1=round(num_th/2)
        x_throat1=geek.linspace(xp_con,L_nozzle_con,num=n1)
        for i in x_throat1:
            th_step_cur=mth.asin((L_nozzle_con-i)/R_u)
            th_step.append(th_step_cur)
            y_cur=R_t+(1-mth.cos(th_step_cur))*R_u
            y_throat1.append(y_cur) # With this part we have defined the coordinates for the first part of the throat

        n2=round(num_th/2)
        x_throat2=geek.linspace(L_nozzle_con,xp,num=n2)
        th_step=[]
        y_throat2=[]
        for i in x_throat2:
            th_step_cur=mth.asin((i-L_nozzle_con)/(Ru_bell*R_t))
            th_step.append(th_step_cur)
            y_cur=R_t+(1-mth.cos(th_step_cur))*R_t*Ru_bell
            y_throat2.append(y_cur)

        n3=num_div
        x_div=geek.linspace(xp,L_tot,num=n3)
        Delta=b**2-4*a*(c-(x_div-L_nozzle_con))
        y_div=geek.zeros(len(x_div))
        for i in range(len(y_div)):
            y_div[i]=(-b+mth.sqrt(Delta[i]))/(2*a) # With this part we have defined the coordinates for the divergent part of the nozzle
        x_2=geek.concatenate((x_throat2,x_div))
        if (n2+n3)>num_cool:
            e_tran_Num=num_cool
            num_th_2=round(len(x_throat2)/len(x_2)*e_tran_Num)
            #A_div_cool=geek.linspace(A_div[0],A_div[-1],e_tran_Num)
            x_th2_cool=geek.linspace(x_throat2[0],x_throat2[-1],num_th_2)
            x_div_cool=geek.linspace(x_div[0],x_div[-1],(e_tran_Num-num_th_2))
            y_th2_cool=[]

            for i in x_th2_cool:
                th_step_cur=mth.asin((i-L_nozzle_con)/(Ru_bell*R_t))
                th_step.append(th_step_cur)
                y_cur=R_t+(1-mth.cos(th_step_cur))*R_t*Ru_bell
                y_th2_cool.append(y_cur);
            Delta_div=b**2-4*a*(c-(x_div_cool-L_nozzle_con))
            y_div_cool=geek.zeros(len(x_div_cool))
            for i in range(len(y_div_cool)):
                y_div_cool[i]=(-b+mth.sqrt(Delta_div[i]))/(2*a)
        else:
            x_div_cool=geek.zeros(len(x_div))
            for i in range(len(x_div_cool)):
                x_div_cool[i]=x_div[i];
            x_th2_cool=geek.zeros(len(x_th2_cool))
            for i in range(len(x_th2_cool)):
                x_th2_cool[i]=x_throat2[i];
            y_th2_cool=geek.zeros(len(y_throat2))
            for i in range(len(y_th2_cool)):
                y_th2_cool[i]=y_throat2[i];
            y_div_cool=geek.zeros(len(y_div))
            for i in range(len(y_div_cool)):
                y_div_cool[i]=y_div[i];
        x_noz=geek.concatenate((x_con, x_throat1, x_throat2, x_div))
        y_noz=geek.concatenate((y_con, y_throat1, y_throat2, y_div));

        x_noz_cool=geek.concatenate((x_con,x_throat1,x_th2_cool,x_div_cool))
        y_noz_cool=geek.concatenate((y_con,y_throat1,y_th2_cool,y_div_cool))
        
    A_noz=mth.pi*y_noz**2
    y2_cool=geek.concatenate((y_th2_cool,y_div_cool))
    A_div_cool=mth.pi*y2_cool**2 # Divergent areas for cooling
    # We are now going to determine the pressure in each section of the nozzle
    
    # Properties in chamber and throat do evaluate properites in the convergent, by taking a linearization from chamber to throat
    Ts=ispObj.get_Temperatures(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state,frozenAtThroat=frozen_state)
    rhos=ispObj.get_Densities(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state,frozenAtThroat=frozen_state)
    cps=ispObj.get_HeatCapacities(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state,frozenAtThroat=frozen_state)
    gs=ispObj.get_Throat_MolWt_gamma(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state)
    gsc=ispObj.get_Chamber_MolWt_gamma(Pc=Pc,MR=MR,eps=eps)

    Transp_c=ispObj.get_Chamber_Transport(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state)
    Transp_t=ispObj.get_Throat_Transport(Pc=Pc,MR=MR,eps=eps,frozen=frozen_state)

    mu_c=Transp_c[1]/10
    mu_t=Transp_t[1]/10

    T_t=Ts[1]
    rho_t=rhos[1]
    rho_c=rhos[0]
    cp_c=cps[0]
    cp_t=cps[1]
    g_t=gs[1]
    g_c=gsc[1]
    R_t=cp_t*(g_t-1)/g_t
    R_c=cp_c*(g_c-1)/g_c

    P_t=rho_t*R_t*T_t # Pressure in the throat

    # For the convergent part we assume a linear decreas in pressure
    x_1=geek.concatenate((x_con,x_throat1))
    a_P_con=(P_t-Pc)/(x_1[-1]-x_1[0])
    P_con=Pc+a_P_con*x_1

    y_2=geek.concatenate((y_throat2,y_div))
    A_div=mth.pi*y_2**2 #Areas sampled at the various locations in the divergent

    P_div=[]
    for i in A_div:
        eps_it=i/At
        dp=ispObj.get_PcOvPe(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state)
        Pe_it=Pc/dp
        P_div.append(Pe_it);
    
    P_noz=geek.concatenate((P_con,P_div)) # Pressure throughout the nozzle

    # We assume that properties change linearly in the convergent

    a_T_con=(T_t-Tc)/(x_1[-1]-x_1[0])
    T_con=Tc+a_T_con*x_1

    T_div=[]
    for i in A_div:
        eps_it=i/At
        T_s=ispObj.get_Temperatures(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state)
        T_it=T_s[2]
        T_div.append(T_it);
    
    T_noz=geek.concatenate((T_con,T_div)) # Temperature throughout the nozzle

    a_R_con=(R_t-R_c)/(x_1[-1]-x_1[0])
    R_con=R_c+a_R_con*x_1 # R of gas in the convergent

    a_g_con=(g_t-g_c)/(x_1[-1]-x_1[0])
    g_con=g_c+a_g_con*x_1 # gamma of gas in the convergent

    a_mu_con=(mu_t-mu_c)/(x_1[-1]-x_1[0])
    mu_con=mu_c+a_mu_con*x_1  # viscosity of gases in the convergent

    a_cp_con=(cp_t-cp_c)/(x_1[-1]-x_1[0])
    cp_con=cp_c+a_cp_con*x_1 # heat capacity in the convergent

    rho_con=P_con/(R_con*T_con)
    u_con=geek.zeros(len(T_con))
    for i in range(len(T_con)):
        u_con[i]=mth.sqrt(T_con[i]*g_con[i]*R_con[i])
    y_1=geek.concatenate((y_con,y_throat1))
    y_2=geek.concatenate((y_throat2,y_div))
    A_con=mth.pi*y_1**2

    v_con=m_p/(rho_con*A_con)
    M_con=v_con/u_con

    Pr_con=4*g_con/(9*g_con-5) # These are all approximations

    r_con=Pr_con**(1/3)

    Tw_ad_con=T_con*(1+r_con*(g_con-1)/2*M_con**2)

    x_2=geek.concatenate((x_throat2,x_div))
    rho_div=[]
    cp_div=[]
    g_div=[]
    v_div=[]
    mu_div=[]
    Pr_div=[]
    k_div=[]
    Mach_div=[]
    g_pr=[]
    
# We loop in the divergent to get the properties in all the sections of the divergent from CEA functions
    for i in A_div_cool:
        eps_it=i/At
        rhos_div=ispObj.get_Densities(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state) # from this the density
        gs_div=ispObj.get_exit_MolWt_gamma(Pc,MR,eps_it)# From this the gamma value
        Mach_n=ispObj.get_MachNumber(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state,frozenAtThroat=frozen_state) # From this the mach number
        Transp=ispObj.get_Exit_Transport(Pc=Pc,MR=MR,eps=eps_it,frozen=frozen_state) # Transport properties (cp, mu, Pr, k)

        rho_div_it=rhos_div[2]
        cp_div_it=Transp[0]
        g_div_it=gs_div[1]
        mu_div_it=Transp[1]/10
        Pr_div_it=Transp[3]
        k_div_it=Transp[2]
        rho_div.append(rho_div_it)
        cp_div.append(cp_div_it)
        g_div.append(g_div_it)
        mu_div.append(mu_div_it)
        Pr_div.append(Pr_div_it)
        Mach_div.append(Mach_n)
        k_div.append(k_div_it);
    
    r_div=geek.zeros(len(Pr_div)) # recovery factor to use in the definition of the adiabatic wall temperature
    Tw_ad_div=geek.zeros(len(Pr_div))
    T_f_div=geek.zeros(len(Pr_div))
    for i in range(len(Pr_div)):
        Pr_act=Pr_div[i]
        r_div[i]=Pr_act**(1/3)
        Tw_ad_div[i]=T_div[i]*(1+r_div[i]*(g_div[i]-1)/2*Mach_div[i]**2)
        T_f_div[i]=0.5*T_w+0.28*T_div[i]+0.22*Tw_ad_div[i]
    Tw_ad_noz=geek.concatenate((Tw_ad_con,Tw_ad_div))

     # Film temperature in the divergent calculated as in Ziebland
    T_f_con=0.5*T_w+0.28*T_con+0.22*Tw_ad_con # Film temperature in convergent

    # Using the Bartz relations we find the convective heat coefficient
    a_b=0.026
    h_c_div=geek.zeros(len(mu_div))
    for i in range(len(mu_div)):
        h_c_div[i]=1.213*a_b*m_p**0.8*mu_div[i]**0.2*cp_div[i]*Pr_div[i]**(-0.6)*(2*y_2[i])**(-1)*(Tc/T_f_div[i])**0.68
    h_c_con=geek.zeros(len(T_f_con))
    for i in range(len(T_f_con)):
        h_c_con[i]=1.213*a_b*m_p**0.8*mu_con[i]**0.2*cp_con[i]*Pr_con[i]**(-0.6)*(2*y_1[i])**(-1)*(Tc/T_f_con[i])**0.68

    h_c_noz=geek.concatenate((h_c_con,h_c_div))

    sig=Material.yieldstress_l
    SF=Default.Safety_factor

    t_noz=SF*P_noz*y_noz/sig #Thickness of the wall in the nozzle
    D_t=R_t*2;
    D_e=mth.sqrt(eps)*D_t

    return t_noz,x_noz,y_noz,Tw_ad_noz,h_c_noz,D_t,D_e,L_nozzle_con,L_nozzle_div,L_tot,x_noz_cool,y_noz_cool,error;




