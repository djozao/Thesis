function convloss = convectionloss_cylinder(Length, Diameter, ApertureDiameter, Temperature, T_amb)
    T_f = (T_amb+Temperature)/2;
    if(T_f >2273.15)
        T_f = 2273.15;
    end
    A_ap = pi*ApertureDiameter^2/4;
    A_cav = pi*Diameter^2/4+pi*Diameter*Length+pi/4*(Diameter^2-ApertureDiameter^2);
    Gr = 9.81*(Temperature-T_amb)*Diameter^3*air.thermalexp(T_f)*air.density(T_f)^2/(air.viscosity(T_f)^2);
    Nu = 7.26*10^(-5)*Gr^(0.3533)*(4+cos(-60*pi/180))^(5.8632)*(A_ap/A_cav)^(1.0266);
    convloss = Nu*A_cav*air.thermalconductivity(T_f)*(Temperature-T_amb)/Diameter; 
end