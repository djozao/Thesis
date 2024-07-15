function [T_insulation, Radloss, Convloss, Condloss] = T_insulation_calculation(L,D,A,t, Tw ,T_amb,emi_ins)
    % This is made considering a cylindrical geometry.
    radloss_frac = outer_radlossfrac_cylinder(L, D, A, emi_ins);
    cond_frac = 2*pi*L/(log((D+2*t)/D));
    %objective function to calculate temperature of the surface
    objective = @(T) abs(cond_frac*(Tw-T)*(0.00665*exp(0.0015*((Tw+T)/2-273.15)))-radloss_frac*(T^4-T_amb^4)-outer_convloss_cylinder(L,D,A,T,T_amb));
    T0 = (T_amb+Tw)/2;
    options =  optimset();
    T_insulation = fminsearch(objective, T0, options);
    Radloss = radloss_frac*(T_insulation^4-T_amb^4);
    Convloss = outer_convloss_cylinder(L,D,A,T_insulation,T_amb);
    Condloss = cond_frac*(Tw-T_insulation)*(0.00665*exp(0.0015*((Tw+T_insulation)/2-273.15)));
end