function  [T_pf, Q_conv, dp, sections] = linear_duct(T_pi, T_w, P_i, A, P, L, m, fluid)
    
    [T_pf1, Q_conv1, dp1] = singlephase_linear_duct(T_pi, T_w, P_i, A, P, L, m, fluid);
    
    if (T_pf1>fluid.boiling(P_i-dp1) && T_pi < fluid.boiling(P_i))
        x = [0,L];
        x_search = (x(1)+x(2))/2;
        while(abs(T_pf1-fluid.boiling(P_i-dp1))>0.01)
            x_search = (x(1)+x(2))/2;
            [T_pf1, Q_conv1, dp1] = singlephase_linear_duct(T_pi, T_w, P_i, A, P, x_search, m, fluid);
            if T_pf1>fluid.boiling(P_i-dp1)
                x(2) = x_search;
            else
                x(1) = x_search;
            end
        end
        [T_pf2, Q_conv2, dp2] = twophase_linear_duct(T_pf1, T_w, P_i-dp1, A, P, L-x_search, m, fluid);
        if(Q_conv2>m*fluid.hvap(P_i-dp1,dp2))
            x = [0,L-x_search];
            while(abs(Q_conv2 - m*fluid.hvap(P_i-dp1,dp2))>0.01)
                x_search2 = (x(1)+x(2))/2;
                [T_pf2, Q_conv2, dp2] = twophase_linear_duct(T_pf1, T_w, P_i-dp1, A, P, x_search2, m, fluid);
                if Q_conv2>m*fluid.hvap(P_i-dp1,dp2)
                    x(2) = x_search2;
                else
                    x(1) = x_search2;
                end
            end
            [T_pf3, Q_conv3, dp3] = singlephase_linear_duct(T_pf2, T_w, P_i-dp1-dp2, A, P, L-x_search-x_search2, m, fluid);
            T_pf = T_pf3;
            Q_conv = Q_conv1+Q_conv2+Q_conv3;
            dp = dp1+dp2+dp3;
            sections= [x_search, x_search2, L-x_search-x_search2];
            disp('Water boiled')
        else
            Q_conv = Q_conv1+Q_conv2;
            dp = dp1+dp2;
            T_pf = T_pf2;
            sections = [x_search,L-x_search,0];
            disp('Water not boiled')

        end

    else
        Q_conv = Q_conv1;
        dp = dp1;
        T_pf= T_pf1;
        sections = [L,0,0];
    end
    
    
end