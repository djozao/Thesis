function T_eq = objective_function(x, CR, P, T_amb, diffuse, Vac)
    T_eq = CalcEqTemp_Geo(CR, P, x(1), x(2), x(3), x(4), x(5), x(6), T_amb, diffuse, Vac);
end