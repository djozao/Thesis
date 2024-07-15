function Power_abs = power_absorbed_function(P, alpha, L, D, A, CR, Abs, diffuse)

inputData = [alpha, Abs, A/D, L/D];
% Convert the numeric array to a table if needed
inputTable = array2table(inputData, 'VariableNames', {'Angle', 'Emi', 'Aperture', 'Length'});

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
        
    else
        disp('choose available cone ratio [0, 0.1, 0.25, 0.5, 0.75, 1')
        return
    end
end