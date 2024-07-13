function [T_pf, Q_prop, dp] = singlephase_linear_duct(T_pi, T_w, P_i, A, P, L, m, fluid)

D_h = 4 * A / P;
T_pf = T_pi;
T_pf_prev = 0;
dp = 0.001;

while (abs(T_pf-T_pf_prev) > 0.1)
    T_pf_prev = T_pf;
    if T_pf-T_pi == 0
        temp_vector = T_pi:1:T_pi;
    else
        temp_vector = T_pi: (T_pf-T_pi)/10:T_pf;
    end
    pressure_vector = P_i:-dp/10:(P_i-dp);
    visc = fluid.viscosity(temp_vector, pressure_vector);
    mean_visc = mean(visc(~isnan(visc)));
    visc(isnan(visc)) = mean_visc;
    visc = mean(visc);
    k = fluid.k( temp_vector, pressure_vector );
    mean_k = mean(k(~isnan(k)));
    k(isnan(k)) = mean_k;
    k = mean(k);
    density = mean(fluid.density( temp_vector, pressure_vector ));
    fluid.cp(temp_vector, pressure_vector);
    cp = mean(fluid.cp(temp_vector, pressure_vector));

    Re_D = 4*m/(visc*P);
    Pr = visc*cp/k;
    if Re_D < 2300
        
        Gz = Re_D*Pr*D_h/L;
        if Gz <= 1000
            Nu_D = 3.657 + 0.2362 * Gz^0.488 * exp(-57.2 / Gz);
        else
            Nu_D = 1.077 * Gz^(1/3) - 0.7;
        end
            h_conv = Nu_D * k / D_h;
            f = 64/Re_D;
            dp = f*L/D_h *0.5*m^2/(A^2*density)*10^-5;
    else
        f = 1 / ((1.82 * log10(Re_D) - 1.64)^2);
        Nu_D = (f / 8) * (Re_D - 1000) * Pr / (1 + 12.7 * sqrt(f / 8) * (Pr^(2/3) - 1)) * (1 + (D_h / L)^(2 / 3));
        h_conv = Nu_D * k / D_h;
        dp = f*L/D_h *0.5*m^2/(A^2*density)*10^-5;
    end
    T_pf = T_pi + (T_w-T_pi) * (1-exp(-h_conv*P*L/(m*cp)));
    %dp = dp + 10^(-5)*((m/A)^2*(1/fluid.density(T_pf,P_i-dp)-1/fluid.density(T_pi,P_i))+9.81*L*1/(visc));
    
     
end
    Q_prop = m*cp*(T_pf-T_pi);
end
