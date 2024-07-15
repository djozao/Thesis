function [T_eq] = CalcEqTemp_FEM_Geo(CR, P, alpha, L, D, A, n_1, n_2,n_3, Abs, emi, t, k, T_amb)
    alpha = alpha*pi/180;
    if A > D
        T_eq = T_amb;
        disp('Aperture cant be bigger than Diameter')
        return;  % Return T_eq = T_amb and exit the function
    end
    o = A/2/tan(alpha);

    if CR == 0
        Gebhart_matrix = Gebhart_3D_Cylinder(L, D, A, Abs, n_1, n_2, n_3);
        Fraction_matrix = fraction_power_matrix_3D_Cylinder(o, L, D, alpha, n_2, n_3);
        Area = zeros(1+n_1+n_2+n_3,1);
        Area(1) = pi*A^2/4;
        for i = 1:n_1
            Area(1+i) = pi/4*( (A+(D-A)/n_1*i)^2 - (A+(D-A)/n_1*(i-1))^2 );
        end
        for i = 1:n_2
            Area(1+n_1+i) = D*L/n_2*pi;
        end
        for i = 1:n_3
            Area(1+n_1+n_2+i) = pi/4*((D/n_3*i)^2-(D/n_3*(i-1))^2);
        end

        Heat = zeros(1+n_1+n_2+n_3,1);
        for i=1:1+n_1+n_2+n_3
            Heat(i) = (1-Abs)*(Gebhart_matrix(1+n_1+1:1+n_1+n_2+n_3,i)')*Fraction_matrix;
            if i> 1+n_1
                Heat(i) = Heat(i) + Abs*Fraction_matrix(i-1-n_1);
            end
        end

        %plot_cylinder(L,D,A,n_1,n_2,n_3,100*Heat(2:size(Heat))./Area(2:size(Area)))
        100*sum(Heat(2:size(Heat)))/sum(Area(2:size(Area)))
        initial_guess = ones(1+n_1 + n_2+ n_3,1 ) * T_amb*4; % All nodes have initial value of T0_initial
        T0_initial = T_amb;
        sum(P*Heat(2:size(Heat)))
        options = optimset('Display', 'off');
        Gebhart_matrix_emi = Gebhart_3D_Cylinder(L, D, A, emi, n_1, n_2, n_3);
        equilibrium_temperature2 = fsolve(@(T) Teq_cylinder_rad(T, n_1, n_2, n_3, emi, L, D, A, Gebhart_matrix_emi, P*Heat(2:size(Heat)), T0_initial), initial_guess, options);
        equilibrium_temperature2 = equilibrium_temperature2(2:size(equilibrium_temperature2,1));
        T_eq = equilibrium_temperature2'*Area(2:length(Area))/sum(Area(2:length(Area)))
        
        %plot_cylinder(L,D,A,n_1,n_2,n_3, equilibrium_temperature2);
        equilibrium_temperature = fsolve(@(T) Teq_cylinder_radcond(T, n_1, n_2, n_3, emi, L, D, A, k, t, Gebhart_matrix_emi, P*Heat(2:size(Heat)), T0_initial), initial_guess, options);
        equilibrium_temperature = equilibrium_temperature(2:size(equilibrium_temperature,1));
        T_eq = equilibrium_temperature'*Area(2:length(Area))/sum(Area(2:length(Area)));
        %plot_cylinder(L,D,A,n_1,n_2,n_3, equilibrium_temperature);

        stef = 5.670374*10^-8;

        T = equilibrium_temperature;
        T2 = equilibrium_temperature2;
        
        %thing = (1/(emi*Gebhart_matrix_emi(2:1+n_1+n_2+n_3,1)'*Area(2:length(Area))))^(1/4);
        %calculates residuals
        res_rad = sum(P*Heat(2:size(Heat)),1) - emi*stef*( (Area(2:1+n_1).*(Gebhart_matrix_emi(2:1+n_1,1)))'*(T2(1:n_1).^4-T_amb^4) +   (Area(1+n_1+1:1+n_1+n_2).*(Gebhart_matrix_emi(1+n_1+1:1+n_1+n_2,1)))'*(T2(n_1+1:n_1+n_2).^4-T_amb^4) + (Area(1+n_1+n_2+1:1+n_1+n_2+n_3).*(Gebhart_matrix_emi(1+n_1+n_2+1:1+n_1+n_2+n_3,1)))'*(T2(n_1+n_2+1:n_1+n_2+n_3).^4-T_amb^4)       )
        res_radcond = sum(P*Heat(2:size(Heat)),1) - emi*stef*( (Area(2:1+n_1).*(Gebhart_matrix_emi(2:1+n_1,1)))'*(T(1:n_1).^4-T_amb^4) +   (Area(1+n_1+1:1+n_1+n_2).*(Gebhart_matrix_emi(1+n_1+1:1+n_1+n_2,1)))'*(T(n_1+1:n_1+n_2).^4-T_amb^4) + (Area(1+n_1+n_2+1:1+n_1+n_2+n_3).*(Gebhart_matrix_emi(1+n_1+n_2+1:1+n_1+n_2+n_3,1)))'*(T(n_1+n_2+1:n_1+n_2+n_3).^4-T_amb^4)       )
 
    elseif CR == 1
        gamma = atan(D/(2*L));
        L_D = L/cos(gamma);
        Gebhart_matrix = Gebhart_3D_Cone(L,D,A,Abs,n_1,n_2);
        Fraction_matrix = fraction_power_3D_cone_matrix(o, L, alpha, gamma, n_2);
        Area = zeros(1+n_1+n_2,1);
        Area(1) = pi*A^2/4;
        
        for i = 1:n_1
            Area(1+i) = pi/4*( (A+(D-A)/n_1*i)^2 - (A+(D-A)/n_1*(i-1))^2 );
        end
        for i = 1:n_2
            Area(1+n_1+i) = (pi*D/2*L_D*(i^2-(i-1)^2)/(n_2^2));
        end

        Heat = zeros(1+n_1+n_2,1);
        for i=1:1+n_1+n_2
            Heat(i) = (1-Abs)*(Gebhart_matrix(1+n_1+1:1+n_1+n_2,i)')*Fraction_matrix;
            if i> 1+n_1
                Heat(i) = Heat(i) + Abs*Fraction_matrix(i-1-n_1);
            end
        end
        sum(P*Heat(2:size(Heat)))
        initial_guess = ones(1+n_1 + n_2,1 ) * T_amb*4; % All nodes have initial value of T0_initial
        T0_initial = T_amb;
        options = optimset('Display', 'off');
        Gebhart_matrix_emi = Gebhart_3D_Cone(L, D, A, emi, n_1, n_2);
        equilibrium_temperature2 = fsolve(@(T) Teq_cone_rad(T, n_1, n_2, emi, L_D, D, A, Gebhart_matrix_emi, P*Heat(2:size(Heat)), T0_initial), initial_guess, options);
        %residuals = Teq_cone_rad(equilibrium_temperature2, n_1, n_2, emi, L_D, D, A, Gebhart_matrix_emi, P * Heat(2:end), T0_initial)
        equilibrium_temperature2 = equilibrium_temperature2(2:size(equilibrium_temperature2,1));
        T_eq = equilibrium_temperature2'*Area(2:length(Area))/sum(Area(2:length(Area)))
        plot_cone(L,D,A,n_1,n_2, equilibrium_temperature2)
        equilibrium_temperature = fsolve(@(T) Teq_cone_radcond(T, n_1, n_2, emi, L_D, D, A, k, t, Gebhart_matrix_emi, P*Heat(2:size(Heat)), T0_initial), initial_guess, options);
        equilibrium_temperature = equilibrium_temperature(2:size(equilibrium_temperature,1));
        T_eq = equilibrium_temperature'*Area(2:length(Area))/sum(Area(2:length(Area)));
        plot_cone(L,D,A,n_1,n_2, equilibrium_temperature)
stef = 5.670374*10^-8;

        T = equilibrium_temperature;
        T2 = equilibrium_temperature2;

        res_rad = sum(P*Heat(2:size(Heat)),1) - emi*stef*( (Area(2:1+n_1).*(Gebhart_matrix_emi(2:1+n_1,1)))'*(T2(1:n_1).^4-T_amb^4) +   (Area(1+n_1+1:1+n_1+n_2).*(Gebhart_matrix_emi(1+n_1+1:1+n_1+n_2,1)))'*(T2(n_1+1:n_1+n_2).^4-T_amb^4) )
        res_radcond = sum(P*Heat(2:size(Heat)),1) - emi*stef*( (Area(2:1+n_1).*(Gebhart_matrix_emi(2:1+n_1,1)))'*(T(1:n_1).^4-T_amb^4) +   (Area(1+n_1+1:1+n_1+n_2).*(Gebhart_matrix_emi(1+n_1+1:1+n_1+n_2,1)))'*(T(n_1+1:n_1+n_2).^4-T_amb^4) )
 

    elseif CR >0
        gamma = atan(D/(2*L*CR));
        L_D = L*CR/cos(gamma);

        Gebhart_matrix = Gebhart_3D_Concyl(L, D, A, Abs, CR, n_1, n_2, n_3);
        Fraction_matrix = fraction_power_matrix_3D_ConCyl(o, L, D, alpha, CR, n_2, n_3);

        Area = zeros(1+n_1+n_2+n_3,1);
        Area(1) = pi*A^2/4;
        for i = 1:n_1
            Area(1+i) = pi/4*( (A+(D-A)/n_1*i)^2 - (A+(D-A)/n_1*(i-1))^2 );
        end
        for i = 1:n_2
            Area(1+n_1+i) = D*L*(1-CR)/n_2*pi;
        end
        for i = 1:n_3
            Area(1+n_1+n_2+i) = (pi*D/2*L_D*(i^2-(i-1)^2)/(n_3^2));
        end

        Heat = zeros(1+n_1+n_2+n_3,1);
        for i=1:1+n_1+n_2+n_3
            Heat(i) = (1-Abs)*(Gebhart_matrix(1+n_1+1:1+n_1+n_2+n_3,i)')*Fraction_matrix;
            if i> 1+n_1
                Heat(i) = Heat(i) + Abs*Fraction_matrix(i-1-n_1);
            end
        end
        sum(P*Heat(2:size(Heat)))
        %plot_concyl(L,D,A,CR,n_1,n_2,n_3,Heat(2:size(Heat)));
        initial_guess = ones(1+n_1 + n_2+n_3,1 ) * T_amb*4; % All nodes have initial value of T0_initial
        T0_initial = T_amb;
        options = optimset('Display', 'off');
        Gebhart_matrix_emi = Gebhart_3D_Concyl(L, D, A, emi, CR, n_1, n_2, n_3);
        equilibrium_temperature2 = fsolve(@(T) Teq_concyl_rad(T, n_1, n_2, n_3, emi, L, D, A, CR, Gebhart_matrix_emi, P*Heat(2:size(Heat)), T0_initial), initial_guess, options);
        equilibrium_temperature2 = equilibrium_temperature2(2:size(equilibrium_temperature2,1));
        T_eq = equilibrium_temperature2'*Area(2:length(Area))/sum(Area(2:length(Area)));
        %plot_concyl(L,D,A,CR,n_1,n_2,n_3, equilibrium_temperature2)


        equilibrium_temperature = fsolve(@(T) Teq_concyl_radcond(T, n_1, n_2, n_3, emi, L, D, A, CR, t, k, Gebhart_matrix_emi, P*Heat(2:size(Heat)), T0_initial), initial_guess, options);
        equilibrium_temperature = equilibrium_temperature(2:size(equilibrium_temperature,1));
        T_eq = equilibrium_temperature'*Area(2:length(Area))/sum(Area(2:length(Area)));
        %plot_concyl(L,D,A,CR,n_1,n_2,n_3, equilibrium_temperature)
        
        stef = 5.670374*10^-8;

        T = equilibrium_temperature;
        T2 = equilibrium_temperature2;

        res_rad = sum(P*Heat(2:size(Heat)),1) - emi*stef*( (Area(2:1+n_1).*(Gebhart_matrix_emi(2:1+n_1,1)))'*(T2(1:n_1).^4-T_amb^4) +   (Area(1+n_1+1:1+n_1+n_2).*(Gebhart_matrix_emi(1+n_1+1:1+n_1+n_2,1)))'*(T2(n_1+1:n_1+n_2).^4-T_amb^4) + (Area(1+n_1+n_2+1:1+n_1+n_2+n_3).*(Gebhart_matrix_emi(1+n_1+n_2+1:1+n_1+n_2+n_3,1)))'*(T2(n_1+n_2+1:n_1+n_2+n_3).^4-T_amb^4)       )
        res_radcond = sum(P*Heat(2:size(Heat)),1) - emi*stef*( (Area(2:1+n_1).*(Gebhart_matrix_emi(2:1+n_1,1)))'*(T(1:n_1).^4-T_amb^4) +   (Area(1+n_1+1:1+n_1+n_2).*(Gebhart_matrix_emi(1+n_1+1:1+n_1+n_2,1)))'*(T(n_1+1:n_1+n_2).^4-T_amb^4) + (Area(1+n_1+n_2+1:1+n_1+n_2+n_3).*(Gebhart_matrix_emi(1+n_1+n_2+1:1+n_1+n_2+n_3,1)))'*(T(n_1+n_2+1:n_1+n_2+n_3).^4-T_amb^4)       )
        %res_radcond_noins = sum(P*Heat(2:size(Heat)),1) -emi*stef*( Area(2:1+n_1+n_2+n_3)'*(T3(1:n_1+n_2+n_3).^4-T_amb^4))-emi*stef*((Area(2:1+n_1+n_2+n_3).*(Gebhart_matrix(2:1+n_1+n_2+n_3,1)))'*(T3(1:n_1+n_2+n_3).^4-293.15^4))       


    else
    
    end




end