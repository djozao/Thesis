function residuals = Teq_cone_rad(T, n_1, n_3, emi, L_D, D, A, B_matrix, Heat, T_0)
    stef = 5.670374*10^-8;
    
    
    % Initialize the residuals vector
    residuals = zeros(n_1 + n_3, 1);
    
    % Calculate residuals for surface 1 nodes
    for j = 1:n_1
        residuals(j) = Heat(j) - emi * pi/4*( (A+(D-A)/n_1*j)^2 - (A+(D-A)/n_1*(j-1))^2 ) * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4));
    end
    
    % Calculate residuals for surface 3 nodes
    for j = 1:n_3
        residuals(n_1 + j) = Heat(n_1 + j) - emi * (pi*D/2*L_D*(j^2-(j-1)^2)/(n_3^2)) * stef * (B_matrix(1 + n_1 + j, :) * (T(1 + n_1 + j)^4 - T.^4));
    end
    residuals = [T_0 - T(1); residuals];
end