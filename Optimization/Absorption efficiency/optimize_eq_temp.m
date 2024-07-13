function optimize_eq_temp()
    % Define constants
    CR = 0; % Example value
    P = 20; % Example value
    T_amb = 50; % Example value
    diffuse = 1; % Example value
    Vac = 1; % Example value

    % Define bounds for the variables [alpha, L, D, A, Abs, emi]
    lb = [5, 0.05, 0.04, 0.039999, 0.5, 0.05];
    ub = [5, 0.1, 0.1, 0.039999, 0.5, 0.05];

    % Define initial guess
    x0 = [5, 0.1, 0.1, 0.39999, 0.5, 0.05];

    % Objective function wrapper
    objective = @(x) -objective_function(x, CR, P, T_amb, diffuse, Vac);

    nonlcon = @(x) nonlinear_constraint(x);

    % Optimization options
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

    % Run the optimizer
    [x_opt, fval] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, options);

    % Display the results
    disp('Optimized variables:');
    disp(['alpha = ', num2str(x_opt(1))]);
    disp(['L = ', num2str(x_opt(2))]);
    disp(['D = ', num2str(x_opt(3))]);
    disp(['A = ', num2str(x_opt(4))]);
    disp(['Abs = ', num2str(x_opt(5))]);
    disp(['emi = ', num2str(x_opt(6))]);
    disp(['Maximum T_eq = ', num2str(-fval)]);
end

