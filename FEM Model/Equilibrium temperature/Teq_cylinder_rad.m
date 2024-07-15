function residuals = Teq_cylinder_rad(T, n_1, n_2, n_3, emi, L, D, A, B_matrix, Heat, T_0)
    stef = 5.670374*10^-8;
    
    % Initialize the residuals vector
    residuals = zeros(n_1 +n_2+ n_3, 1);
    
    % Calculate residuals for surface 1 nodes
    for j = 1:n_1
        residuals(j) = Heat(j) - emi * ((A+(D-A)/n_1*j)^2-(A+(D-A)/n_1*(j-1))^2)*pi /4 * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4));
    end

    % Calculate residuals for surface 1 nodes
    for j = 1:n_2
        residuals(n_1+j) = Heat(n_1+j) - emi * L / n_2 * pi * D * stef * (B_matrix(1 + n_1 + j, :) * (T(1+n_1+j)^4 - T.^4));
    end
    
    % Calculate residuals for surface 3 nodes
    for j = 1:n_3
        residuals(n_1+n_2 + j) = Heat(n_1 + n_2 + j) - emi * ((D*j/n_3)^2-(D*(j-1)/n_3)^2)* pi / 4 * stef * (B_matrix(1 + n_1 + n_2+ j, :) * (T(1 + n_1 + n_2+ j)^4 - T.^4));
    end
    residuals = [T_0 - T(1); residuals];
end