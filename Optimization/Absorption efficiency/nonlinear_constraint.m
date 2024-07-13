% Nonlinear constraint function (D >= A)
function [c, ceq] = nonlinear_constraint(x)
    c = x(4) - x(3); % D >= A
    ceq = []; % No equality constraint
end