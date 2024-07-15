function rad_loss = radiationlossfrac_cone(Length, Diameter, ApertureDiameter, Emissivity) 
    L = Length/(Diameter/2);
    R = ApertureDiameter/Diameter;
    F11 = 1- 1/(sqrt(1+L^2));
    F12 = (1-R^2)/(sqrt(1+L^2));
    F13 = R^2 /(sqrt(1+L^2));
    F21 = 1;
    F22 = 0;
    F23 = 0;

    syms B13 B23;
    eqs= [B13 == F13 + (1-Emissivity)*(F11*B13+F12*B23),B23 == F23 + (1-Emissivity)*(F21*B13+F22*B23)];
    Bsolve = solve(eqs, [B13, B23]);
    rad_loss = double(Emissivity*5.670374419*10^(-8)*(Bsolve.B13*pi*Diameter/2*sqrt(Diameter^2/4+Length^2)+Bsolve.B23*pi*(Diameter^2-ApertureDiameter^2)/4)); 
end 