function [Ts, TC, ms, mD, Tra_s_vec, Tra_r_vec,Es_vec, Cap_flux_vec, f_E] = The_Model_f_E_vec(F, T_R, q_R, P,NDVI,steps_per_day)
                               % radiaiton Temp humidity precip              
                             
    %%% PHYSICAL CONSTANTS %%%

    N = length(F);      % Number of soil moisture values we simulate
    rho_a = 1.25;       % density of air [kg/m^3]
    c_p = 1003;         % specific heat of air [J/kg/K]    
    c_s = 10000;         % specific heat of dry soil [J/kg/K]
    rho_s = 1000;       % density of dry soil [kg/m^3]
    L = 2257000;        % Latent enthalpy of vaporization [J/kg]
    T_freeze = 273.15;  % Freezing temperature [K]
    P_s = 900;          % surface pressure [hPa]
    rho_l = 1000;       % density of water [kg/m^3]    

    %%% TUNABLE PARAMETERS (or ones I'm not confident about) %%%

      % c_leaf =  1000;     % specific heat of the canopy surface [J/kg/K]  % ORIGINAL
      % rho_leaf = 300;     % density of leaf canopy [kg/m^3]


        c_leaf =  10000;     % specific heat of the canopy surface [J/kg/K]
        rho_leaf = 300;     % density of leaf canopy [kg/m^3]


                 g_s = 1/2000;       % dry surface conductance [m/s] ORIGINAL
                 g_C = 1/800;        % vegetation conductance from canopy [m/s]
       

         % g_s = 1/1000;       % dry surface conductance [m/s] used for 02. run
         % g_C = 1/400;        % vegetation conductance from canopy [m/s]

        % g_s = 1/4000;       % dry surface conductance [m/s] 
        % g_C = 1/1600;        % vegetation conductance from canopy [m/s]


      % alpha_s = 0.75;     % Surface albedo
      % alpha_C = 0.9;      % Plant albedo

       alpha_s = 1.0;     % Surface albedo
       alpha_C = 1.0;      % Plant albedo

       % alpha_s = 0.9;     % Surface albedo
       % alpha_C = 1;      % Plant albedo


       %nu_H = 3;          % Dry surface sensitivity
       %nu_C = 3;  


       % MJB comment .. both were 20 originally. But runs were with 2.
       nu_H = 2;          % Dry surface sensitivity
       nu_C = 2;      

    tau_cap = pi * 1e7; % Capillary timescale (1 year)



     % minmax scaling to 0-1
     % f_E = (NDVI - (-0.2)) ./ (1 - (-0.2)) ;  
      % f_E = (NDVI - min(NDVI)) ./ (max(NDVI) - min(NDVI) ) ; %+ randn(1,9131) .* 0.07 ;


     % f_E = 0.3 .* sin(linspace(0,N/365*pi,N)) + 0.45 + randn(1,N) * 0.07 ; 


     % f_E = 0.3 .* sin(linspace(0,N/365*pi,N)) + 0.45 + randn(1,N) * 0.07 ; 

     % MJB comment: Lets assume max realistic range of NDVI is roughly
     % 0.8-0.9. We want to scale NDVI dynamics relative to that range, but
     % still allow for pixels that have NDVI dynamics +- 0. With minmax
     % scaling they get set to 0 or 1 if their mean is too high low.
     % Therefore do minmax scaling of the range rather than everything.
     % Then translate to f_e with min of 0 


       f_E = (NDVI - (0.0)) ./ (0.8 - (0.0)) ;   
       f_E(f_E < 0) = 0 ; 
       f_E(f_E > 1) = 1 ; 





    r = 0.7;            % Root fraction in surface layer

    %%% GEOMETRY %%%
    % MJB comments could decrease h_s .. should lead to bigger impact of
    % f_E

    theta_max = 0.6;    % soil pore space [-]
    % h_s = 0.1;          % meter
    h_s = 0.1;          % meter
    h_d = 1;            % meter
    hC = 1;             % meter

    %%% Time step and other initializations %%%

    sec_per_day = 86400;    % seconds per day
    dt = 86400 / steps_per_day;    % time increment (10 chunks per day)
    i = 1;

    Ts = zeros(1, N); % temperature soil
    ms = zeros(1, N); % moisture surface
    mD = zeros(1, N); % moisture deep
    TC = zeros(1, N); % temperature canopy
    Tra_s_vec = zeros(1, N); % transpiration shallow
    Tra_r_vec = zeros(1, N); % transpiration root
    Cap_flux_vec = zeros(1, N); % capilariy flux to shallow     
    Es_vec = zeros(1, N); % evaporation shallow   


    % initialize with non zeors
    Ts(1) = T_R(1);
    ms(1) = 0.3;
    mD(1) = 0.3;
    TC(1) = T_R(1);




    %%% Running the model %%%

    while i <  N

        %%% Sensible Heat Fluxes %%%

        H_s = nu_H * (Ts(i) - T_R(i));
        H_C = nu_C * (TC(i) - T_R(i));

        %%% Surface Evaporation %%%

        qss = calc_q_s(Ts(i) - T_freeze, P_s);          % Saturation Specific Humidity
        qs_def = qss - q_R(i);                          % Specific Humidity Gradient
        if qs_def > 0
            Es = rho_a * g_s * ms(i) * qs_def;          % Evapotranspiration
        else
            Es = 0;
        end

        %%% Transpiration %%%

        qsC = calc_q_s(TC(i) - T_freeze, P_s);
        qC_def = qsC - q_R(i);

        if qC_def > 0
            Tra_s = rho_a * g_C * ms(i) * qC_def * r;          % Canopy transpiration from the surface layer
            Tra_r = rho_a * g_C * mD(i) * qC_def * (1 - r);    % Same from the root layer
        else
            Tra_s = 0;
            Tra_r = 0;
        end

        %%% Capillary flux %%%

        if ms(i) < mD(i)
            Cap_flux = rho_l * h_s * (mD(i) - ms(i)) / tau_cap;    % Capillary action acts only when the surface is drier than depth

        else
            Cap_flux = 0;
        end
            % Cap_flux = 0;

        %%% ENERGY BUDGETS %%%

        dTs_dt = (alpha_s * (1 - f_E(:,i)) * F(i) - H_s - L * Es) / (c_s * rho_s * h_s);
        dTC_dt = (alpha_C * f_E(:,i) * F(i) - H_C - L * (Tra_s + Tra_r)) / (c_leaf * rho_leaf * hC);
        % dTs_dt = (alpha_s * (1 - f_E) * F(i) - H_s - L * Es) / (c_s * rho_s * h_s);
        % dTC_dt = (alpha_C * f_E * F(i) - H_C - L * (Tra_s + Tra_r)) / (c_leaf * rho_leaf * hC);


        %%% MOISTURE BUDGETS %%%

        dms_dt = (P(i) - Es - Tra_s + Cap_flux) / (rho_l * h_s * theta_max);
        dmd_dt = -(Tra_r + Cap_flux) / (rho_l * h_d * theta_max);

        Tra_s_vec(i + 1) = Tra_s ; 
        Tra_r_vec(i + 1) = Tra_r ; 
        Es_vec(i + 1) =  Es ; 
        Cap_flux_vec(i + 1) = Cap_flux ;        

        %%% INTEGRATING %%%

        Ts(i + 1) = Ts(i) + dTs_dt * dt;
        TC(i + 1) = TC(i) + dTC_dt * dt;
        ms(i + 1) = ms(i) + dms_dt * dt;
        mD(i + 1) = mD(i) + dmd_dt * dt;



        %%% Overflow into the deep layer %%%
        % while there is overflow into deeper layer which is kinda like
        % drainage, there is no real mechanims for higher loss if both
        % shallow and deep layers are full. So "drainage" can only happen
        % if deep is not full, which is rare, normally deep is full.

        if ms(i + 1) > mD(i + 1)
            ro = (ms(i + 1) - mD(i + 1)) * h_s / h_d;
            ms(i + 1) = mD(i + 1);
            mD(i + 1) = mD(i + 1) + ro;
        end
        if ms(i + 1) < 0
            ms(i + 1) = 0;
        end
        if mD(i + 1) > 1
            mD(i + 1) = 1;
        end
        if mD(i + 1) < 0
            mD(i + 1) = 0;
        end
        i = i + 1;


    end

end





%% plots for diagnosis


         % figure
         % plot(Ts(1:i-1),'b.-') ; hold on ; plot(TC(1:i-1),'r.-')  ;
         % legend('Tsoil','Tcanopy')  
         % 
         % figure
         % plot(ms(1:i-1),'b.-') ; hold on ; plot(mD(1:i-1),'r.-')  ;
         % legend('SM surf','SM deep')  








