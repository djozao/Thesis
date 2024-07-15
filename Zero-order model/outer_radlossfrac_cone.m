function radlossfrac = outer_radlossfrac_cone(Length, Diameter, ApertureDiameter, Emissivity)
    radlossfrac = Emissivity*5.670374419*10^(-8)*(pi/4*(Diameter^2-ApertureDiameter^2)+pi*Diameter/2*sqrt(Length^2+Diameter^2/4));
end