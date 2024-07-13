function T_eq = CalcEqTemp_Geo(CR, P, alpha, L, D, A, Abs, emi, T_amb, diffuse, Vac)
    % Determine power absorbed depending on surface conditions

    if A > D
        T_eq = T_amb;
        Power_abs = 0;
        disp('Aperture cant be bigger than Diameter')
        return;  % Return T_eq = T_amb and exit the function
    end
    stefan = 5.670374*10^-8;
    inputData = [alpha, Abs, A/D, L/D];
    % Convert the numeric array to a table if needed
    inputTable = array2table(inputData, 'VariableNames', {'Angle', 'Emi', 'Aperture', 'Length'});
    
    % You'll need to change the location depend on your file organization
    if CR == 0
        if diffuse == 0
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\spec_CC00.mat");
            Power_abs = P*spec_CC00.predictFcn([alpha,Abs,A/D,L/D])/100;
        elseif diffuse == 1
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\dif_CC00.mat")
            Power_abs = P*dif_CC00.predictFcn(inputTable)/100;
            
        else
            disp('Choose specular or diffuse reflection (0, 1)')
        end
        
        Gebhart_matrix = Gebhart_3D_Cylinder(L,D,A,emi, 1, 1, 1);
        Area = [pi*(D^2-A^2)/4, pi*D*L, pi*D^2/4];
        Rad_frac = stefan*emi*(Area*Gebhart_matrix(2:4,1));
    elseif CR == 0.1
        if diffuse == 0
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\spec_CC01.mat");
            Power_abs = P*spec_CC01.predictFcn([alpha,Abs,A/D,L/D])/100;
        elseif diffuse == 1
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\dif_CC01.mat");
            Power_abs = P*dif_CC01.predictFcn(inputTable)/100;
        else
            disp('Choose specular or diffuse reflection (0, 1)')
        end
         
        Gebhart_matrix = Gebhart_3D_Concyl(L,D,A,emi,CR, 1, 1, 1);
        Area = [pi*(D^2-A^2)/4, pi*D*L*(1-CR), pi*D/2*sqrt((D/2)^2+(CR*L)^2  )];
        Rad_frac = stefan*emi*(Area*Gebhart_matrix(2:4,1));
    elseif CR == 0.25
        if diffuse == 0
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\spec_CC025.mat");
            Power_abs = P*spec_CC025.predictFcn([alpha,Abs,A/D,L/D])/100;
        elseif diffuse == 1
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\dif_CC025.mat");
            Power_abs = P*dif_CC025.predictFcn(inputTable)/100;
        else
            disp('Choose specular or diffuse reflection (0, 1)')
        end
        Gebhart_matrix = Gebhart_3D_Concyl(L,D,A,emi,0.25, 1, 1, 1);
        Area = [pi*(D^2-A^2)/4, pi*D*L*(1-CR), pi*D/2*sqrt((D/2)^2+(CR*L)^2  )];
        Rad_frac = stefan*emi*(Area*Gebhart_matrix(2:4,1));
        
    elseif CR == 0.5
        if diffuse == 0
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\spec_CC05.mat");
            Power_abs = P*spec_CC05.predictFcn([alpha,Abs,A/D,L/D])/100;
        elseif diffuse == 1
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\dif_CC05.mat");
            Power_abs = P*dif_CC05.predictFcn(inputTable)/100;
        else
            disp('Choose specular or diffuse reflection (0, 1)')
        end
        Gebhart_matrix = Gebhart_3D_Concyl(L,D,A,emi,0.5, 1, 1, 1);
        Area = [pi*(D^2-A^2)/4, pi*D*L*(1-CR), pi*D/2*sqrt((D/2)^2+(CR*L)^2  )];
        Rad_frac = stefan*emi*(Area*Gebhart_matrix(2:4,1));
        
    elseif CR == 0.75
        if diffuse == 0
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\spec_CC075.mat");
            Power_abs = P*spec_CC075.predictFcn([alpha,Abs,A/D,L/D])/100;
        elseif diffuse == 1
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\dif_CC075.mat");
            Power_abs = P*dif_CC075.predictFcn(inputTable)/100;
        else
            disp('Choose specular or diffuse reflection (0, 1)')
        end
        Gebhart_matrix = Gebhart_3D_Concyl(L,D,A,emi,0.75, 1, 1, 1);
        Area = [pi*(D^2-A^2)/4, pi*D*L*(1-CR), pi*D/2*sqrt((D/2)^2+(CR*L)^2  )];
        Rad_frac = stefan*emi*(Area*Gebhart_matrix(2:4,1));

    elseif CR == 1
        if diffuse == 0
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\spec_CC1.mat");
            Power_abs = P*spec_CC1.predictFcn([alpha,Abs,A/D,L/D])/100;
        elseif diffuse == 1
            load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\dif_CC1.mat");
            Power_abs = P*dif_CC1.predictFcn(inputTable)/100;
        else
            disp('Choose specular or diffuse reflection (0, 1)')
        end
        Gebhart_matrix = Gebhart_3D_Cone(L,D,A,emi, 1, 1);
        Area = [pi*(D^2-A^2)/4, pi*D/2*sqrt((D/2)^2+(L)^2  )];
        Rad_frac = stefan*emi*(Area*Gebhart_matrix(2:3,1));
    else
        disp('choose available cone ratio [0, 0.1, 0.25, 0.5, 0.75, 1')
        return
    end
    
    if Power_abs < P*Abs 
        Power_abs = P*Abs;
    elseif Power_abs > P
        Power_abs = P;
    end
    %add radation/convection outer loss in the equations if wanted
    if Vac == 1
        T_eq = (Power_abs/Rad_frac+T_amb^4)^(1/4);
    elseif Vac == 0
        %define convection loss depending on geometry in the function
        objective = @(T) abs(Power_abs-(T^4-T_amb^4)*Rad_frac-convectionloss_cylinder(L,D,A,T,T_amb));
        T0 = T_amb;
        options =  optimset();
        T_eq = fminsearch(objective, T0, options);
        if (T_eq+T_amb)/2 > 2273.15
            fprintf("Air has surpassed its temperature limit; Results might be unrealistic")
        end
    end
    %Rad_frac*(T_eq^4-T_amb^4)
    %Power_abs
end