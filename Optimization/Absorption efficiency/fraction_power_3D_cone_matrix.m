function fraction_power_ring = fraction_power_3D_cone_matrix(o, L, alpha, gamma, n_3)

L_D = L/cos(gamma);
h_max = tan(alpha)*(o+L)/(1+tan(alpha)/tan(gamma));
fraction_power_ring = zeros(n_3,1);

j = 0;
for i = 1:n_3
    if sin(gamma)*L_D/n_3*i < h_max
        fraction_power_ring(i,1) = (atan((sin(gamma)*L_D/n_3*i)/(o+L-(sin(gamma)*L_D/n_3*i)/tan(gamma)))^2-atan((sin(gamma)*L_D/n_3*(i-1))/(o+L-(sin(gamma)*L_D/n_3*(i-1))/tan(gamma)))^2)/(alpha^2);
    elseif j == 0
        fraction_power_ring(i,1) = (atan((h_max)/(o+L-(h_max)/tan(gamma)))^2-atan((sin(gamma)*L_D/n_3*(i-1))/(o+L-(sin(gamma)*L_D/n_3*(i-1))/tan(gamma)))^2)/(alpha^2);
        j = 1;
    else
        fraction_power_ring(i,1) = 0;
    end
end