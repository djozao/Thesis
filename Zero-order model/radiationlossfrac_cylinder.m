function rad_loss = radiationlossfrac_cylinder(Length, Diameter, ApertureDiameter, Emissivity) 
    R1 = Diameter/(2*Length);
    R2 = ApertureDiameter/(2*Length);
    x = 1 + (1+R2^2)/R1^2;
    F11 = 0;
    F14 = 1/2*(x-sqrt(x^2-4*(ApertureDiameter/Diameter)^2));
    F13 = 1/2*(2+1/R1-sqrt((2+1/R1)^2-4))-F14;
    F12 = 1-F13-F14;
    F21 = R1/2*F12;
    F22 = 1-2*F21;
    F24 = 1/(Diameter*Length)*(ApertureDiameter^2/4-Diameter^2/4*F14);
    F23 = 1-F21-F22-F24;
    if(ApertureDiameter~=Diameter)
        F31 = F13*1/(1-(ApertureDiameter/Diameter)^2);
        F32 = 1- F31;
    else 
        F31 = 0;
        F32 = 0;
    end
    F33 = 0;
    F34 = 0;
               
    syms B14 B24 B34;
    eqs= [B14 == F14 + (1-Emissivity)*(F11*B14+F12*B24+F13*B34),B24 == F24 + (1-Emissivity)*(F21*B14+F22*B24+F23*B34), B34 == F34 + (1-Emissivity)*(F31*B14+F32*B24+F33*B34)];
    Bsolve = solve(eqs, [B14, B24, B34]);
    rad_loss = double(Emissivity*5.670374419*10^(-8)*(Bsolve.B14*pi*Diameter^2/4+Bsolve.B24*Diameter*Length*pi+Bsolve.B34*pi*(Diameter^2-ApertureDiameter^2)/4));
end 