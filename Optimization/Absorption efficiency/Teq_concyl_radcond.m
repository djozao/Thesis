function residuals = Teq_concyl_radcond(T, n_1, n_2, n_3, emi, L, D, A, CR, t, k, B_matrix, Heat, T_0)
    stef = 5.670374*10^-8;
    gamma = atan(D/(2*L*CR));
    L_D = L*CR/cos(gamma);
    % Initialize the residuals vector
    residuals = zeros(n_1 +n_2+ n_3, 1);
    
    % Calculate residuals for surface 1 nodes
    for j = 1:n_1
        if n_1 == 1
            residuals(j) = Heat(j) -...
            emi * ((A+(D-A)/n_1*j)^2-(A+(D-A)/n_1*(j-1))^2)*pi /4 * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4))...
            + 2*pi*k*t/(log(  D/(A+(D-A)/2)  ))*(T(1+j+n_2)-T(1+j))...
            + k*pi/4*((D+t)^2-D^2)*( T(1+j+n_2)-T(1+j))/(L*(1-CR)/2/n_2);
           
        elseif j == 1

            residuals(j) = Heat(j) -...
            emi * ((A+(D-A)/n_1*j)^2-(A+(D-A)/n_1*(j-1))^2)*pi /4 * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4))...
            + 2*pi*k*t/(log(  (A+(D-A)/n_1/2+(D-A)/n_1*(j)) /(A+(D-A)/n_1/2+(D-A)/n_1*(j-1))  ))*(T(1+j+1)-T(1+j));

        elseif j<n_1
            residuals(j) = Heat(j) -...
            emi * ((A+(D-A)/n_1*j)^2-(A+(D-A)/n_1*(j-1))^2)*pi /4 * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4))...
            + 2*pi*k*t/(log(  (A+(D-A)/n_1/2+(D-A)/n_1*(j)) /(A+(D-A)/n_1/2+(D-A)/n_1*(j-1))  ))*(T(1+j+1)-T(1+j))...
            + 2*pi*k*t/(log(  (A+(D-A)/n_1/2+(D-A)/n_1*(j-1)) /(A+(D-A)/n_1/2+(D-A)/n_1*(j-2)) ))*(T(1+j-1)-T(1+j));
        else
            
            %residuals(j) = Heat(j) - emi * ((A+(D-A)/n_1*j)^2-(A+(D-A)/n_1*(j-1))^2)*pi /4 * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4));
            residuals(j) = Heat(j) -...
            emi * ((A+(D-A)/n_1*j)^2-(A+(D-A)/n_1*(j-1))^2)*pi /4 * stef * (B_matrix(1 + j, :) * (T(1+j)^4 - T.^4))...
            + 2*pi*k*t/(log(  D/(D-(D-A)/n_1/2)  ))*(T(1+j+n_2)-T(1+j))...
            + k*pi/4*((D+t)^2-D^2)*( T(1+j+n_2)-T(1+j))/(L*(1-CR)/2/n_2)...
            + 2*pi*k*t/(log((D-(D-A)/n_1/2)/(D-(D-A)/n_1/2-(D-A)/n_1)   ))*(T(1+j-1)-T(1+j));

        end
    end

    % Calculate residuals for surface 1 nodes
    for j = 1:n_2
        if n_2 == 1
            residuals(n_1+j) = Heat(n_1+j) - emi * L*(1-CR) / n_2 * pi * D * stef * (B_matrix(1 + n_1 + j, :) * (T(1+n_1+j)^4 - T.^4))...
                + 2*pi*k*t/(log( D/(D-(D-A)/n_1/2) ))*(T(1+n_1)-T(1+n_1+n_2))...
                + k*pi/4*((D+t)^2-D^2)*(T(1+n_1)-T(1+n_1+n_2))/(L*(1-CR)/2/n_2)...
                + (T(1+n_1+n_2+n_3)-T(1+j+n_1))* 2*pi*k*t*tan(gamma)/(log(D/2   /  (D/2-sin(gamma)*L_D/2/n_3)))...
                + k*pi/4*((D+t)^2-D^2)*(T(1+n_1+n_2+n_3)-T(1+n_1+n_2))/(L*(1-CR)/2/n_2);
         elseif j==1
             residuals(n_1+j) = Heat(n_1+j) - emi * L*(1-CR) / n_2 * pi * D * stef * (B_matrix(1 + n_1 + j, :) * (T(1+n_1+j)^4 - T.^4))...
                 + k*pi/4*((D+t)^2-D^2)*(T(1+n_1+j+1)-T(1+n_1+j))/(L*(1-CR)/n_2)...
                 + (T(1+n_1+n_2+n_3)-T(1+j+n_1))* 2*pi*k*t*tan(gamma)/(log(D/2   /  (D/2-sin(gamma)*L_D/2/n_3)))...
                 + k*pi/4*((D+t)^2-D^2)*(T(1+n_1+n_2+n_3)-T(1+n_1+j))/(L*(1-CR)/2/n_2);
         elseif j<n_2
             residuals(n_1+j) = Heat(n_1+j) - emi * L*(1-CR) / n_2 * pi * D * stef * (B_matrix(1 + n_1 + j, :) * (T(1+n_1+j)^4 - T.^4))...
                 + k*pi/4*((D+t)^2-D^2)*(T(1+n_1+j+1)-T(1+n_1+j))/(L*(1-CR)/n_2) ...
                 + k*pi/4*((D+t)^2-D^2)*(T(1+n_1+j-1)-T(1+n_1+j))/(L*(1-CR)/n_2) ...
                 ;
         else
             residuals(n_1+j) = Heat(n_1+j) - emi * L*(1-CR) / n_2 * pi * D * stef * (B_matrix(1 + n_1 + j, :) * (T(1+n_1+j)^4 - T.^4))...
                 + 2*pi*k*t/(log( D/(D-(D-A)/n_1/2) ))*(T(1+n_1)-T(1+n_1+j))...
                 + k*pi/4*((D+t)^2-D^2)*(T(1+n_1)-T(1+n_1+j))/(L*(1-CR)/2/n_2) +...
                 + k*pi/4*((D+t)^2-D^2)*(T(1+n_1+j-1)-T(1+n_1+j))/(L*(1-CR)/n_2);
        end
    end
    

    
    for j = 1:n_3
        if n_3 == 1
            residuals(n_1+n_2 + j) = Heat(n_1+n_2 + j) - emi * (pi*D/2*L_D*(j^2-(j-1)^2)/(n_3^2)) * stef * (B_matrix(1 + n_1 + n_2+ j, :) * (T(1 + n_1+ n_2+ + j)^4 - T.^4))...
                                + k*pi/4*((D+t)^2-D^2)*(T(1+n_1+1)-T(1+n_1+n_2+j))/(L*(1-CR)/2/n_2)...
                                + (T(1+n_1+1)-T(1+n_2+j+n_1))* 2*pi*k*t*tan(gamma)/(log(D/2   /  (D/2-sin(gamma)*L_D/2/n_3)));

        elseif j == 1
            residuals(n_1+n_2 + j) = Heat(n_1+n_2 + j) - emi * (pi*D/2*L_D*(j^2-(j-1)^2)/(n_3^2)) * stef * (B_matrix(1 + n_1 +  n_2+ j, :) * (T(1 + n_1+ n_2 + j)^4 - T.^4))...
                                + (T(1+n_1+n_2+j+1)- T(1+n_1+n_2+j)) * 2*pi*k*t*tan(gamma) / ( log( (sin(gamma)*L_D/2/n_3 + sin(gamma)*L_D/n_3)/(sin(gamma)*L_D/2/n_3) ) );

        elseif j == n_3

            residuals(n_1+n_2 + j) = Heat(n_1+n_2 + j) - emi * (pi*D/2*L_D*(j^2-(j-1)^2)/(n_3^2)) * stef * (B_matrix(1 + n_1 + n_2 +j, :) * (T(1 + n_1+ n_2 + j)^4 - T.^4))...
                                + k*pi/4*((D+t)^2-D^2)*(T(1+n_1+1)-T(1+n_1+n_2+j))/(L*(1-CR)/2/n_2)...
                                + (T(1+n_1+1)-T(1+j+n_1+n_2))* 2*pi*k*t*tan(gamma)/(log(D/2   /  (D/2-sin(gamma)*L_D/2/n_3)))...
                                + (T(1+n_1+n_2+j-1)-T(1+n_2+j+n_1))* 2*pi*k*t*tan(gamma)/(log( (D/2-sin(gamma)*L_D/2/n_3) / (D/2-sin(gamma)*L_D/2/n_3-sin(gamma)*L_D/n_3)));
        else

            residuals(n_1 +n_2+ j) = Heat(n_1+n_2 + j) - emi * (pi*D/2*L_D*(j^2-(j-1)^2)/(n_3^2)) * stef * (B_matrix(1 + n_1 + n_2 + j, :) * (T(1 + n_1+ n_2 + j)^4 - T.^4))...
                + (T(1+n_1+n_2+j-1)-T(1+j+n_2+n_1))* 2*pi*k*t*tan(gamma)/(log( (sin(gamma)*L_D/2/n_3 + sin(gamma)*L_D/n_3*(j-1)) / (sin(gamma)*L_D/2/n_3 + sin(gamma)*L_D/n_3*(j-2) )  ))...
                + (T(1+n_1+n_2+j+1)-T(1+j+n_2+n_1))* 2*pi*k*t*tan(gamma)/(log( (sin(gamma)*L_D/2/n_3 + sin(gamma)*L_D/n_3*(j)) / (sin(gamma)*L_D/2/n_3 + sin(gamma)*L_D/n_3*(j-1) )  ));

        end
    end
    residuals = [T_0 - T(1); residuals];
end