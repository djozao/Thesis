function convloss = outer_convloss_cone(L, D, A, T, T_amb)
    T_f = (T+T_amb)/2;
    slant = sqrt(L^2+D^2/4);
    Gr = 9.81*air.thermalexp(T_f)*air.density(T_f)^2*(T-T_amb)*slant^3/(air.viscosity(T_f)^2);
    Nu = 0.65*((Gr*cos(2*arctan(D/(2*L))))^0.25+1.44/(D/(2*L)));
    h = Nu*air.thermalconductivity(T_f)/slant;
    convloss = h*(T-T_amb)*(pi/4*(D^2-A^2)+pi*D/2*sqrt(L^2+D^2/4));
end