function percentage = concyl_fem_absorbed(alpha, L, D, A, CR, emi, n_1, n_2, n_3)
alpha = alpha*pi/180;
o = A/2/tan(alpha);
Gebhart_matrix = Gebhart_3D_Concyl(L, D, A, emi, CR, n_1, n_2,n_3);
Fraction_matrix = fraction_power_matrix_3D_ConCyl(o, L, D, alpha, CR, n_2, n_3);
Heat = zeros(1+n_1+n_2+n_3,1);
%Calculates the heat flux
for i=1:1+n_1+n_2+n_3
    Heat(i) = (1-emi)*(Gebhart_matrix(1+n_1+1:1+n_1+n_2+n_3,i)')*Fraction_matrix;
    if i> 1+n_1
        Heat(i) = Heat(i) + emi*Fraction_matrix(i-1-n_1);
    end
end
percentage= 100*(1-Heat(1));

end