%Define constants
P_i = 5;
L = 0.1;
number_pipes = 2;
T_p = 293.15;
T_w = 293.15;
T_f = 700; % bulk temperature
fluid = nitrogen(T_f:1:T_f+1, P_i-1:1:P_i+1);
% Define initial guesses and bounds for optimization variables
min_height = 0.001;
max_height = 0.003;
min_side = 0.0005;
max_side = 0.005;

H_guess = (min_height+max_height)/2;  % Initial guess for h
W_guess = (min_side+max_side)/2;  % Initial guess for w
m_guess = 251*10^-6;  % Initial guess for m
initial_guess = [H_guess, W_guess, m_guess];

lb = [ min_height , min_side , m_guess  ];  % Lower bounds for optimization variables
ub = [ max_height , max_side , m_guess  ];  % Upper bounds for optimization variables

% Define anonymous function for optimization
fun = @(x) -linear_duct_efficiency(P_i, x(1), x(2), L, x(3), number_pipes, fluid, T_f);

% Perform optimization using patternsearch
options = optimoptions('patternsearch','Display','off');
[x_opt, fval] = patternsearch(fun, initial_guess, [], [], [], [], lb, ub, [], options);

% Display results
disp('Optimal Parameters:');
disp(x_opt);
disp('Optimal Efficiency:');
disp(-fval);

