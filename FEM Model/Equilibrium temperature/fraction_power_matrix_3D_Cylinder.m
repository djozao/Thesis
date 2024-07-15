function fraction_power_ring = fraction_power_matrix_3D_Cylinder(o, L, D, alpha, n_2, n_3)

fraction_power_ring = zeros(n_2+n_3,1);
alpha_sep = atan((D/2)/(o+L));

if alpha_sep > alpha
    h_max = tan(alpha)*(o+L);
    fraction_power_ring(1:n_2) = 0;
    j=0;
    for i = 1:n_3
        if D/n_3/2*i < h_max
            fraction_power_ring(n_2+i) =( (1-cos(atan((D/n_3/2*i)/(o+L)))-(1-cos(atan((D/n_3/2*(i-1))/(o+L))))) )/(1-cos(alpha));
        elseif j == 0
            fraction_power_ring(n_2+i) = (  1-cos(atan((h_max)/(o+L)))-  (1-cos(atan((D/n_3/2*(i-1))/(o+L))))  )/(1-cos(alpha) );
            j = 1;
        else
            fraction_power_ring(n_2+i) = 0;
        end
    end
else
    for i = 1:n_3
        fraction_power_ring(n_2+i) = (   1-cos(atan((D/n_3/2*i)/(o+L)))- (1-cos(atan((D/n_3/2*(i-1))/(o+L)))) )/(1-cos(alpha));        
    end
    x_max = D/(2*tan(alpha))-o;
    j = 0;
    for i = 1:n_2
        if L/n_2*i < L-x_max
            fraction_power_ring(i) = (   1-cos(atan((D/2)/(o+L-L/n_2*i)))-(1-cos(atan((D/2)/(o+L-L/n_2*(i-1))))))/(1-cos(alpha));
        elseif j == 0
            fraction_power_ring(i) = (1-cos((alpha))-(1-cos(atan((D/2)/(o+L-L/n_2*(i-1))))))/(1-cos(alpha));
            j = 1;
        else
            fraction_power_ring(i) = 0;
        end
    end
end