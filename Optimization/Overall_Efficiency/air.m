classdef air
    methods(Static)

        function [density] = density(Temperature)

            x = Temperature-273.15;
            Temp = [0,5,10,15,20,25,30,40,50,60,80,100,125,150,175,200,225,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,2000]; % in degrees Celsius
            density_table = [1.292,1.268,1.246,1.225,1.204,1.184,1.164,1.127,1.093,1.060,1,0.9467,0.8868,0.8338,0.7868,0.7451,0.7078,0.6168,0.5238,0.4567,0.4043,0.3526,0.3289,0.3009,0.2773,0.2571,0.2395,0.2242,0.2109,0.1992,0.1885,0.1553]; % in kg/m^3
            if x>2000
                density = 28.9647/(8.31446261815324*Temperature)*101325/1000;
            else
                density = interp1(Temp,density_table,x,"linear"); % your equation or look-up table
            end 
        end

        function [thermalexp] = thermalexp(Temperature)

            x = Temperature-273.15;
            Temp = [0,5,10,15,20,25,30,40,50,60,80,100,125,150,175,200,225,300,400,500,600,700,800,900,1000,1100]; % in degrees Celsius
            thermalexp_table = [3.69,3.62,3.56,3.50,3.43,3.38,3.32,3.21,3.12,3.02,2.85,2.70,2.51,2.33,2.22,2.10,2.01,1.76,1.52,1.32,1.16,1.03,0.94,0.86,0.80,0.75]; % in kg/m^3
            if (x<=1100)
                thermalexp = interp1(Temp,thermalexp_table,x,"linear")*10^(-3); % your equation or look-up table
            else
                thermalexp= 1/(Temperature);
            end
        end

        function [viscosity] = viscosity(Temperature)

            x = Temperature-273.15;
            Temp = [0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000]; % in degrees Celsius
            viscosity_table = [1.729,1.754,1.778,1.802,1.825,1.849,1.872,1.895,1.918,1.941,1.963,2.008,2.052,2.096,2.139,2.181,2.264,2.345,2.420,2.504,2.577,2.760,2.934,3.101,3.261,3.415,3.563,3.846,4.111,4.362,4.600,4.826,5.817,6.630];
            viscosity = interp1(Temp,viscosity_table,x,"linear")*10^(-5); % your equation or look-up table
        end
        function [thermalcond] = thermalconductivity(Temperature)

            x = Temperature-273.15;
            Temp = [0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1500,2000]; % in degrees Celsius
            thermalconductivity_table = [0.02364,0.02401,0.02439,0.02476,0.02514,0.02551,0.02588,0.02625,0.02662,0.02699,0.02735,0.02808,0.02881,0.02953,0.03024,0.03095,0.03235,0.03374,0.03511,0.03646,0.03779,0.04104,0.04418,0.04721,0.05015,0.05298,0.05572,0.06093,0.06581,0.07037,0.07465,0.07868,0.09599,0.11113];
            thermalcond = interp1(Temp,thermalconductivity_table,x,"linear"); % your equation or look-up table
        end
        function [heatcapacity] = cp(Temperature)
            Temp = [60, 78.79, 81.61, 100, 120,140, 160,180,200,220,300,320,340,360,380,400,500,600,700,800,900,1100,1500,1900];
            cp_table = [1.901,1.933, 1.089,1.040,1.022,1.014,1.011,1.008,1.007,1.006,1.006,1.007,1.009,1.010,1.012,1.014,1.030,1.051,1.075,1.099,1.121,1.159,1.210,1.241];
            heatcapacity = interp1(Temp, cp_table,Temperature,"linear")*10^3;
        end


    end 
end