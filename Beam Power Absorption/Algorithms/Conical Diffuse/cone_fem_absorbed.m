function percentage = cone_fem_absorbed(alpha, L, D, A, emi, n_1, n_2)
alpha = alpha*pi/180;
o = A/2/tan(alpha);
gamma = atan(D/(2*L));

Gebhart_matrix = Gebhart_3D_Cone(L, D, A, emi, n_1, n_2);
Fraction_matrix = fraction_power_3D_cone_matrix(o, L, alpha,gamma, n_2);
Heat = zeros(1+n_1+n_2,1);
%Calculates the heat flux
for i=1:1+n_1+n_2
    Heat(i) = (1-emi)*(Gebhart_matrix(1+n_1+1:1+n_1+n_2,i)')*Fraction_matrix;
    if i> 1+n_1
        Heat(i) = Heat(i) + emi*Fraction_matrix(i-1-n_1);
    end
end
percentage= 100*(1-Heat(1));

end