function fraction_power_ring = fraction_power_matrix_3D_ConCyl(o, L, D, alpha, CC, n_2, n_3)



fraction_power_ring = zeros(n_2+n_3,1);
L_C = L*CC;
alpha_sep = atan((D/2)/(o+L*(1-CC)));
gamma = atan(D/(2*L_C));
L_D = L_C/cos(gamma);
L_cyl = L*(1-CC);

if alpha_sep > alpha
    h_max = tan(alpha)*(o+L)/(1+tan(alpha)/tan(gamma));
    fraction_power_ring(1:n_2) = 0;
    j=0;
    for i = 1:n_3
        if sin(gamma)*L_D/n_3*i < h_max
            fraction_power_ring(n_2+i,1) = ( 1-cos(atan((sin(gamma)*L_D/n_3*i)/(o+L-(sin(gamma)*L_D/n_3*i)/tan(gamma))))-(1-cos(atan((sin(gamma)*L_D/n_3*(i-1))/(o+L-(sin(gamma)*L_D/n_3*(i-1))/tan(gamma))))))/(1-cos(alpha));

        elseif j == 0
            fraction_power_ring(n_2+i,1) = ( 1-cos(atan((h_max)/(o+L-(h_max)/tan(gamma))))-(1-cos(atan((sin(gamma)*L_D/n_3*(i-1))/(o+L-(sin(gamma)*L_D/n_3*(i-1))/tan(gamma))))))/(1-cos(alpha));
            j = 1;
        else
            fraction_power_ring(n_2+i,1) = 0;
        end
    end
else
    for i = 1:n_3
        fraction_power_ring(n_2+i,1) = (atan((sin(gamma)*L_D/n_3*i)/(o+L-(sin(gamma)*L_D/n_3*i)/tan(gamma)))^2-atan((sin(gamma)*L_D/n_3*(i-1))/(o+L-(sin(gamma)*L_D/n_3*(i-1))/tan(gamma)))^2)/(alpha^2);
    end
    x_max = D/(2*tan(alpha))-o;
    j = 0;
    for i = 1:n_2
        if L_cyl/n_2*i < L_cyl-x_max
            fraction_power_ring(i,1) = (1-cos(atan((D/2)/(o+L_cyl-L_cyl/n_2*i)))-(1-cos(atan((D/2)/(o+L_cyl-L_cyl/n_2*(i-1))))))/(1-cos(alpha));
        elseif j == 0
            fraction_power_ring(i,1) = (1-cos(alpha)-(1-cos(atan((D/2)/(o+L_cyl-L_cyl/n_2*(i-1))))))/(1-cos(alpha));
            j = 1;
        else
            fraction_power_ring(i,1) = 0;
        end
    end
end