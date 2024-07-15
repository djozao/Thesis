function [T_pf, Q_conv, dp] = singlephase_spiral_duct(T_pi, T_w, P_i, A, P, L, D_c, m, fluid)

D_h = 4 * A / P;
T_pf = T_w;
T_pf_prev = 0;
dp = 0.5;
Re_D_critical = 20000*(D_h/D_c);

while (abs(T_pf-T_pf_prev) > 0.001)
    
    T_pf_prev = T_pf;
    temp_vector = T_pi: (T_pf-T_pi)/10000:T_pf;
    pressure_vector = P_i:-dp/10000:(P_i-dp);
    visc = mean(fluid.viscosity(temp_vector, pressure_vector ));
    k = mean(fluid.k( temp_vector, pressure_vector ));
    density = mean(fluid.density( temp_vector, pressure_vector ));
    cp = mean(fluid.cp(temp_vector, pressure_vector ));

    Re_D = 4*m/(visc*P);
    Pr = visc*cp/k;
    if Re_D < Re_D_critical
            Nu_D = 0.913*(Re_D*(D_h/D_c)^0.5)^0.476*Pr^0.2;
            h_conv = Nu_D * k/ D_h;
            f = 64/Re_D * (1/ (1-(1-(11.6/(Re_D*(D_h/D_c)^0.5)^0.45)^(1/0.45))));
            dp = f*L/D_h *0.5*m^2/(A^2*density)*10^-5;
    else
        f = 0.304*Re_D^(-0.25)+0.029*(D_h/D_c)^(0.5);
        Nu_D = 0.023*Re_D^0.85*Pr^0.4*(D_h/D_c)^0.1;
        h_conv = Nu_D * k / D_h;
        dp = f*L/D_h *0.5*m^2/(A^2*density)*10^-5;
    end
    T_pf = T_pi + (T_w-T_pi) * (1-exp(-h_conv*P*L/(m*cp)));
end
dp
    Q_conv = m*cp*(T_pf-T_pi);
end