function residuals = Teq_cone_radcond(T, n_1, n_2, emi, L_D, D, A, k, t, B_matrix, Heat, T_0)
    stef = 5.670374*10^-8;
    
    gamma = asin(D/2/L_D);
    % Initialize the residuals vector
    residuals = zeros(n_1 + n_2, 1);
    
    for j = 1:n_1
        if n_1 == 1
            residuals(j) = Heat(j) - emi * pi/4*( (A+(D-A)/n_1*j)^2 - (A+(D-A)/n_1*(j-1))^2 ) * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4))...
            + 2*pi*k*t/(log(  D/(A+(D-A)/2)  ))*(T(1+j+n_2)-T(1+j))...
            + (T(1+j+n_2)-T(1+j))* 2*pi*k*t*tan(gamma)/(log( (D/2)   /  (D/2-sin(gamma)*L_D/2/n_2)));

        elseif j == 1
            residuals(j) = Heat(j) - emi * pi/4*( (A+(D-A)/n_1*j)^2 - (A+(D-A)/n_1*(j-1))^2 ) * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4))...
                        + 2*pi*k*t/(log(  (A+(D-A)/n_1/2+(D-A)/n_1*(j)) /(A+(D-A)/n_1/2+(D-A)/n_1*(j-1))  ))*(T(1+j+1)-T(1+j));
        
        elseif j == n_1
            residuals(j) = Heat(j) - emi * pi/4*( (A+(D-A)/n_1*j)^2 - (A+(D-A)/n_1*(j-1))^2 ) * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4))...
                        + 2*pi*k*t/(log(  D/(D-(D-A)/n_1/2)  ))*(T(1+j+n_2)-T(1+j))...  
                        + (T(1+j+n_2)-T(1+j))* 2*pi*k*t*tan(gamma)/(log( (D/2)   /  (D/2-sin(gamma)*L_D/2/n_2)))...
                        + 2*pi*k*t/(log((D-(D-A)/n_1/2)/(D-(D-A)/n_1/2-(D-A)/n_1)   ))*(T(1+j-1)-T(1+j));
            
        else

        residuals(j) = Heat(j) - emi * pi/4*( (A+(D-A)/n_1*j)^2 - (A+(D-A)/n_1*(j-1))^2 ) * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4))...
            + 2*pi*k*t/(log(  (A+(D-A)/n_1/2+(D-A)/n_1*(j)) /(A+(D-A)/n_1/2+(D-A)/n_1*(j-1))  ))*(T(1+j+1)-T(1+j))...
            + 2*pi*k*t/(log(  (A+(D-A)/n_1/2+(D-A)/n_1*(j-1)) /(A+(D-A)/n_1/2+(D-A)/n_1*(j-2)) ))*(T(1+j-1)-T(1+j));
        end
    end


    
    % Calculate residuals for surface 3 nodes
    for j = 1:n_2
        if n_2 == 1
            residuals(n_1 + j) = Heat(n_1 + j) - emi * (pi*D/2*L_D*(j^2-(j-1)^2)/(n_2^2)) * stef * (B_matrix(1 + n_1 + j, :) * (T(1 + n_1 + j)^4 - T.^4))...
                                + 2*pi*k*t/(log(  D/(D-(D-A)/n_1/2)  ))*(T(1+n_1)-T(1+n_1+j))...
                                + (T(1+n_1)-T(1+j+n_1))* 2*pi*k*t*tan(gamma)/(log((D/2)   /  (D/2-sin(gamma)*L_D/2/n_2)));

        elseif j == 1
            residuals(n_1 + j) = Heat(n_1 + j) - emi * (pi*D/2*L_D*(j^2-(j-1)^2)/(n_2^2)) * stef * (B_matrix(1 + n_1 + j, :) * (T(1 + n_1 + j)^4 - T.^4))...
                                + (T(1+n_1+j+1)- T(1+n_1+j)) * 2*pi*k*t*tan(gamma) / ( log( (sin(gamma)*L_D/2/n_2 + sin(gamma)*L_D/n_2)/(sin(gamma)*L_D/2/n_2) ) );

        elseif j == n_2

            residuals(n_1 + j) = Heat(n_1 + j) - emi * (pi*D/2*L_D*(j^2-(j-1)^2)/(n_2^2)) * stef * (B_matrix(1 + n_1 + j, :) * (T(1 + n_1 + j)^4 - T.^4))...
                                + 2*pi*k*t/(log(  D/(D-(D-A)/n_1/2)  ))*(T(1+n_1)-T(1+n_1+j))...
                                + (T(1+n_1)-T(1+j+n_1))* 2*pi*k*t*tan(gamma)/(log(D/2   /  (D/2-sin(gamma)*L_D/2/n_2)))...
                                + (T(1+n_1+j-1)-T(1+j+n_1))* 2*pi*k*t*tan(gamma)/(log( (D/2-sin(gamma)*L_D/2/n_2) / (D/2-sin(gamma)*L_D/2/n_2-sin(gamma)*L_D/n_2)));
        else


            residuals(n_1 + j) = Heat(n_1 + j) - emi * (pi*D/2*L_D*(j^2-(j-1)^2)/(n_2^2)) * stef * (B_matrix(1 + n_1 + j, :) * (T(1 + n_1 + j)^4 - T.^4))...
                + (T(1+n_1+j-1)-T(1+j+n_1))* 2*pi*k*t*tan(gamma)/(log( (sin(gamma)*L_D/2/n_2 + sin(gamma)*L_D/n_2*(j-1)) / (sin(gamma)*L_D/2/n_2 + sin(gamma)*L_D/n_2*(j-2) )  ))...
                + (T(1+n_1+j+1)-T(1+j+n_1))* 2*pi*k*t*tan(gamma)/(log( (sin(gamma)*L_D/2/n_2 + sin(gamma)*L_D/n_2*(j)) / (sin(gamma)*L_D/2/n_2 + sin(gamma)*L_D/n_2*(j-1) )  ));

        end
    end
    residuals = [T_0 - T(1); residuals];
end