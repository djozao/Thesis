function convloss = outer_convloss_cylinder(L, D, A, T, T_amb)
    T_f = (T+T_amb)/2;
    Gr = 9.81*air.thermalexp(T_f)*air.density(T_f)^2*(T-T_amb)*D^3/(air.viscosity(T_f)^2);
    Pr = air.viscosity(T_f)*air.cp(T_f)/air.thermalconductivity(T_f);
    Ra = Gr*Pr;
    Nu = sqrt(0.6 + 0.387*(Ra/((1+(0.559/Pr)^(9/16))^(16/9)))^(1/6));
    h = Nu*air.thermalconductivity(T_f)/D;
    convloss = h*(T-T_amb)*(2*pi*D^2/4-pi*A^2/4+pi*L*D);
end