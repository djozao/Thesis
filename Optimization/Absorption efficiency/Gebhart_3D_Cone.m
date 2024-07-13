function Gebhart_matrix = Gebhart_3D_Cone(L, D, A, emi, n_1, n_3)

gamma = atan(D/(2*L));
L_D = L/cos(gamma);
vf_matrix = zeros(1+n_1+n_3,1+n_1+n_3);

vf_circlecircle = @(r1, r2, h) 1/2*((1 + (1+(r2/h)^2)/((r1/h)^2)) - sqrt((1 + (1+(r2/h)^2)/((r1/h)^2))^2-4*(r2/r1)^2));

% F0_3i
for i = 1:n_3
    if i == n_3
        vf_matrix(1, 1+n_1+i) = 1- vf_circlecircle(A/2, sin(gamma)*L_D/n_3*(i-1), L-cos(gamma)*L_D/n_3*(i-1));
    elseif i == 1
        vf_matrix(1, 1+n_1+i) = vf_circlecircle(A/2, sin(gamma)*L_D/n_3, L-cos(gamma)*L_D/n_3);
        %for j = 1:n_3
        %    if j == 1
        %        vf_matrix(1+n_1+i, 1+n_1+i) = 1 - pi*(sin(gamma)*L_D/n_3)/(pi*L_D/n_3);  
        %    else
        %        vf_matrix(1+n_1+i, 1+n_1+j) = - sum(vf_matrix(1+n_1+i, 1+n_1+2:1+n_1+j-1))+ pi*(sin(gamma)*L_D/n_3)/(pi*L_D/n_3) - vf_circlecircle(sin(gamma)*L_D/n_3*j, sin(gamma)*L_D/n_3, L_D*cos(gamma)/n_3*(j-1)) * (pi*(sin(gamma)*L_D/n_3*j)^2)/(pi*D/2*L_D*(i^2-(i-1)^2)/(n_3^2));
        %        sum(vf_matrix(1+n_1+i, 1+n_1+2:1+n_1+j-1))
        %    end
        %end
    else
        vf_matrix(1, 1+n_1+i) = vf_circlecircle(A/2, sin(gamma)*L_D/n_3*i, L-cos(gamma)*L_D/n_3*i)-vf_circlecircle(A/2, sin(gamma)*L_D/n_3*(i-1), L-cos(gamma)*L_D/n_3*(i-1));
        
    end
        vf_matrix(1+n_1+i, 1) = vf_matrix(1, 1+n_1+i)* (pi*A^2/4) / (pi*D/2*L_D*(i^2-(i-1)^2)/(n_3^2));

end

% F1_3i
for i = 1:n_1
    for j = 1:n_3
        if j == n_3
            vf_matrix(1+i, 1+n_1+j) = (1+ (A/2+((D-A)/(2*n_1)*(i-1)))^2/((A/2+((D-A)/(2*n_1)*(i)))^2 -(A/2+((D-A)/(2*n_1)*(i-1)))^2))* ( 1-vf_circlecircle( (A/2+((D-A)/(2*n_1)*(i))), sin(gamma)*L_D*(j-1)/n_3, L-cos(gamma)*L_D/n_3*(j-1))) - (A/2+((D-A)/(2*n_1)*(i-1)))^2/((A/2+((D-A)/(2*n_1)*(i)))^2 -(A/2+((D-A)/(2*n_1)*(i-1)))^2)*(1 -vf_circlecircle( (A/2+((D-A)/(2*n_1)*(i-1))), sin(gamma)*L_D*(j-1)/n_3, L-cos(gamma)*L_D/n_3*(j-1))); 
        else
            vf_matrix(1+i, 1+n_1+j) = (1+ (A/2+((D-A)/(2*n_1)*(i-1)))^2/((A/2+((D-A)/(2*n_1)*(i)))^2 -(A/2+((D-A)/(2*n_1)*(i-1)))^2))* ( vf_circlecircle( (A/2+((D-A)/(2*n_1)*(i))), sin(gamma)*L_D*j/n_3, L-cos(gamma)*L_D/n_3*j) -vf_circlecircle( (A/2+((D-A)/(2*n_1)*(i))), sin(gamma)*L_D*(j-1)/n_3, L-cos(gamma)*L_D/n_3*(j-1))) - (A/2+((D-A)/(2*n_1)*(i-1)))^2/((A/2+((D-A)/(2*n_1)*(i)))^2 -(A/2+((D-A)/(2*n_1)*(i-1)))^2)*( vf_circlecircle( (A/2+((D-A)/(2*n_1)*(i-1))), sin(gamma)*L_D*j/n_3, L-cos(gamma)*L_D/n_3*j) -vf_circlecircle( (A/2+((D-A)/(2*n_1)*(i-1))), sin(gamma)*L_D*(j-1)/n_3, L-cos(gamma)*L_D/n_3*(j-1))); 
        end
            vf_matrix(1+n_1+j, 1+i) = vf_matrix(1+i, 1+n_1+j) * ((A/2+((D-A)/(2*n_1)*(i)))^2 -(A/2+((D-A)/(2*n_1)*(i-1)))^2)/((D/2*L_D*(j^2-(j-1)^2)/(n_3^2)));
    end
end


% F3i_3j
for i = 1:n_3
    for j = i:n_3
        if i == 1
            if j == 1
                vf_matrix(1+n_1+i, 1+n_1+i) = 1 - pi*(sin(gamma)*L_D/n_3)/(pi*L_D/n_3);  
            else
                vf_matrix(1+n_1+i, 1+n_1+j) = - sum(vf_matrix(1+n_1+i, 1+n_1+2:1+n_1+j-1))+ pi*(sin(gamma)*L_D/n_3)/(pi*L_D/n_3) - vf_circlecircle(sin(gamma)*L_D/n_3*j, sin(gamma)*L_D/n_3, L_D*cos(gamma)/n_3*(j-1)) * (pi*(sin(gamma)*L_D/n_3*j)^2)/(pi*D/2*L_D*(i^2-(i-1)^2)/(n_3^2));
                vf_matrix(1+n_1+j, 1+n_1+i) = vf_matrix(1+n_1+i, 1+n_1+j)* (i^2-(i-1)^2)/(j^2-(j-1)^2);
            %sum(vf_matrix(1+n_1+i, 1+n_1+2:1+n_1+j-1))
            end
        else
            if i == j
                %vf_matrix(1+n_1+i, 1+n_1+j) = ( sin(gamma)*(L_D*i/n_3)^2*(1-sin(gamma)) -((D/2*L_D*(j^2-(j-1)^2)/(n_3^2)))*sum(vf_matrix(1+n_1+i, 1+n_1+1:1+n_1+i-1)) - ((D/2*L_D*((j-1)^2)/(n_3^2)))*(vf_matrix(1+n_1+i-1,1+n_1+i-1) +vf_matrix(1+n_1+i-1, 1+n_1+i))) / ((D/2*L_D*(j^2-(j-1)^2)/(n_3^2)));
                vf_matrix(1+n_1+i, 1+n_1+j) = (sin(gamma)*(L_D*i/n_3)^2*(1-sin(gamma)) -((D/2*L_D*(j^2-(j-1)^2)/(n_3^2)))*sum(vf_matrix(1+n_1+i, 1+n_1+1:1+n_1+i-1)) - ((D/2*L_D*((j-1)^2)/(n_3^2)))*( 1-sin(gamma) + sin(gamma) - vf_circlecircle(sin(gamma)*L_D/n_3*i, sin(gamma)*L_D/n_3*(i-1), L_D*cos(gamma)/n_3) * (pi*(sin(gamma)*L_D/n_3*(j))^2)/(pi*D/2*L_D*((i-1)^2)/(n_3^2))      )  ) / ((D/2*L_D*(j^2-(j-1)^2)/(n_3^2)));
            else
                 %               vf_matrix(1+n_1+i, 1+n_1+j) = ( ((D/2*L_D*(j^2)/(n_3^2)))*aux_self -((D/2*L_D*(j^2-(j-1)^2)/(n_3^2)))*sum(vf_matrix(1+n_1+i, 1+n_1+1:1+n_1+i-1))- ((D/2*L_D*((j-1)^2)/(n_3^2)))*(vf_matrix(1+n_1+i-1,1+n_1+i-1) +vf_matrix(1+n_1+i-1, 1+n_1+i))  ) / ((D/2*L_D*(j^2-(j-1)^2)/(n_3^2)));
                A_ij = pi*sin(gamma)*(L_D/n_3*i)^2;
                A_j = pi*sin(gamma)*(L_D/n_3*(i-1))^2;
                if j == i+1
                    F_ij_b = sin(gamma);
                else
                    F_ij_b = vf_circlecircle(sin(gamma)*L_D/n_3*(j-1), sin(gamma)*L_D/n_3*i, L_D/n_3*cos(gamma)*(j-i-1))* (pi*(sin(gamma)*L_D/n_3*(j-1))^2)/(pi*L_D/n_3*i*sin(gamma)*L_D/n_3*i);
                end
                F_ij_a = vf_circlecircle(sin(gamma)*L_D/n_3*(j), sin(gamma)*L_D/n_3*i, L_D/n_3*cos(gamma)*(j-i))* (pi*(sin(gamma)*L_D/n_3*(j))^2)/(pi*L_D/n_3*i*sin(gamma)*L_D/n_3*i);
                F_j_b = vf_circlecircle(sin(gamma)*L_D/n_3*(j-1), sin(gamma)*L_D/n_3*(i-1), L_D/n_3*cos(gamma)*(j-i))* (pi*(sin(gamma)*L_D/n_3*(j-1))^2)/(pi*L_D/n_3*(i-1)*sin(gamma)*L_D/n_3*(i-1));
                F_j_a = vf_circlecircle(sin(gamma)*L_D/n_3*(j), sin(gamma)*L_D/n_3*(i-1), L_D/n_3*cos(gamma)*(j-i+1))* (pi*(sin(gamma)*L_D/n_3*(j))^2)/(pi*L_D/n_3*(i-1)*sin(gamma)*L_D/n_3*(i-1));
                A_i = A_ij-A_j;
                vf_matrix(1+n_1+i, 1+n_1+j) = (A_ij*(F_ij_b-F_ij_a)-A_j*(F_j_b-F_j_a))/A_i;
                vf_matrix(1+n_1+j, 1+n_1+i) = vf_matrix(1+n_1+i, 1+n_1+j)*A_i/(pi*D/2*L_D*(j^2-(j-1)^2)/(n_3^2));
            end
           
        end

    end
end
%b = -vf_matrix(1:1+n_1+n_3,1);
A_matrix = vf_matrix*(1-emi);
A_matrix(:,1) = A_matrix(:,1)*0;
A_matrix = A_matrix - eye(1+n_3+n_1);
%X = linsolve(A_matrix,b);

B = -vf_matrix;
B(:,2:1+n_1+n_3) = B(:,2:1+n_1+n_3)*emi;
Gebhart_matrix = linsolve(A_matrix, B);