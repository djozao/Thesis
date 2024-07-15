function [T_pf, Q_conv, dp] = twophase_linear_duct(T_pi, T_w, P_i, A, P, L, m, fluid)

    D_h = 4 * A / P;
    dp = 0.01;
    dp_prev = 0;
    while (abs(dp-dp_prev)>1*10^-4)
        dp_prev = dp;

        [viscL, viscG] = fluid.Bviscosity(P_i, dp);
        visc = (viscL+viscG)/2;
        [kL, kG] = fluid.Bk(P_i,dp);
        k = (kL+kG)/2;
        [densityL, densityG] = fluid.Bdensity(P_i,dp);
        density = (densityL+densityG)/2;
        [cpL, cpG] = fluid.Bcp(P_i,dp);
        cp = (cpL+cpG)/2;

        Re_D = 4*m/(visc*P);
        Pr = visc*cp/k;
        if Re_D < 2300
            
                f = 64/Re_D;
                dp = f*L/D_h *0.5*m^2/(A^2*densityL)*10^-5;
        else
            f = 1 / ((1.82 * log10(Re_D) - 1.64)^2);
            
            dp = f*L/D_h *0.5*m^2/(A^2*densityL)*10^-5;
        end
        h_conv = 0.023*k/D_h*(Re_D*density/densityG)^0.8*Pr^0.4;
        X = (densityG/densityL)^0.5*(viscL/viscG)^0.1;
        dp = dp*(1+20/X+1/(X^2));
        dp_mom = (m/A)^2*(1/densityG-1/densityL)*10^-5;
        dp = dp + dp_mom;
    end
    T_pf = fluid.boiling(P_i-dp);
    Q_conv = h_conv*P*L*(T_w-T_pi);
    
end