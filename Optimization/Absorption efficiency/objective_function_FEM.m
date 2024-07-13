function T_eq = objective_function_FEM(x, CR, P, n_1, n_2, n_3, k, T_amb)
    T_eq = CalcEqTemp_FEM_Geo(CR, P, x(1), x(2), x(3), x(4), n_1, n_2, n_3, x(5), x(6), x(7), k, T_amb);
end