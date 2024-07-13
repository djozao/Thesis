function optimize_eq_temp_FEM()
    % Define constants
    CR = 0.5; % Example value
    P = 20; % Example value
    T_amb = 50; % Example value
    k = 237;
    t = 0.004;
    n_1 = 10;
    n_2 = 10;
    n_3 = 10;
    % Define bounds for the variables [alpha, L, D, A, Abs, emi]
    %lb = [5, 1, 1, 1, 0.1, 0.1];
    %ub = [20, 100, 100, 1, 0.4, 0.7];

    lb = [5, 0.05, 0.04, 0.039999, 0.15, 0.05, t];
    ub = [5, 0.1, 0.1, 0.039999, 0.15, 0.05, t];

    % Define initial guess
    x0 = [5, 0.1, 0.1, 0.39999, 0.15, 0.15, t];

    % Objective function wrapper
    objective = @(x) -objective_function_FEM(x, CR, P, n_1, n_2, n_3, k, T_amb);

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