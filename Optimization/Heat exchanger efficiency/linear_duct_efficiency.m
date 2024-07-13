function efficiency = linear_duct_efficiency(P_i, height, width, L, m, number_pipes, fluid, T_f)
    P = (height+width)*2;
    A = height*width;
    D_h = 4 * A / P;
    dp = 0.1;
    dp_prev = 0;
    m = m/number_pipes;
    while (abs(dp-dp_prev) > 0.0001)
        dp_prev = dp;
        
        temp_vector = T_f:1:T_f;
        
       
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
                Nu_D = 3.657/(tanh(2.264*Gz^(-1/3)+1.7*Gz^(-2/3)))+0.0499*Gz*tanh(1/Gz);
                h_conv = Nu_D * k / D_h;
                f = 64/Re_D;
                dp = f*L/D_h *0.5*m^2/(A^2*density)*10^-5;
        else
            f = 1 / ((1.82 * log10(Re_D) - 1.64)^2);
            Nu_D = (f / 8) * (Re_D - 1000) * Pr / (1 + 12.7 * sqrt(f / 8) * (Pr^(2/3) - 1)) * (1 + (D_h / L)^(2 / 3));
            h_conv = Nu_D * k / D_h;
            dp = f*L/D_h *0.5*m^2/(A^2*density)*10^-5;
            disp('turb')
        end
    end
    efficiency = 1- exp(-h_conv*P*L/(m*cp));
end