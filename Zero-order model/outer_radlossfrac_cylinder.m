function radlossfrac = outer_radlossfrac_cylinder(Length, Diameter, ApertureDiameter, Emissivity)
    radlossfrac = Emissivity*5.670374419*10^(-8)*(2*pi*Diameter^2/4+pi*Diameter*Length-pi*ApertureDiameter^2/4);
end