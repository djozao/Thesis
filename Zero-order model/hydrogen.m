classdef hydrogen
    properties (SetAccess = public)
        data
        data2
    end
    methods
        function obj = hydrogen(Temperature, Pressure)
            if nargin > 0
                 [obj.data, obj.data2] = nistdata('H2',Temperature,Pressure);
            end
        end
        function [viscosity] = viscosity(obj, Temperature, Pressure) %Pa s
            viscosity = interp2(obj.data.P,obj.data.T,obj.data.mu,Pressure*10^5,Temperature);
            checkNaN(viscosity, 'viscosity', Temperature, Pressure);
        end
        function [cp] = cp(obj, Temperature, Pressure) % J/kg/K
            cp = interp2(obj.data.P,obj.data.T,obj.data.Cp,Pressure*10^5,Temperature)/obj.data.Mw;
            checkNaN(cp, 'Cp', Temperature, Pressure);
        end
        function [k] = k(obj, Temperature, Pressure)    % W/m/K
            k = interp2(obj.data.P,obj.data.T,obj.data.k,Pressure*10^5,Temperature);
            checkNaN(k, 'Thermal Conductivity', Temperature, Pressure);
        end
        function density = density(obj,Temperature, Pressure) % kg/m^3
            density = interp2(obj.data.P,obj.data.T,obj.data.Rho,Pressure*10^5,Temperature)*obj.data.Mw;
            checkNaN(density, 'Density', Temperature, Pressure);
        end
        function enthalpy = enthalpy(obj,Temperature, Pressure) % kg/m^3
            enthalpy = interp2(obj.data.P,obj.data.T,obj.data.H,Pressure*10^5,Temperature)/obj.data.Mw;
            checkNaN(enthalpy, 'Enthalpy', Temperature, Pressure);
        end
        function T_boiling = boiling(obj, Pressure)
            T_boiling = interp1(obj.data2.P(1,:),obj.data2.T(1,:), Pressure*10^5);
            checkNaN(T_boiling, 'Boiling Temperature',[] , Pressure);
        end
        function hvap = hvap(obj, Pressure, dp)
            hvap = (interp1(obj.data2.P(2,:),obj.data2.H(2,:), Pressure*10^5-dp)- interp1(obj.data2.P(1,:),obj.data2.H(1,:), Pressure*10^5))/obj.data.Mw;
            checkNaN(hvap, 'Hvaporization', Pressure, Pressure-dp);
        end
        function [viscliquid, viscgas] = Bviscosity(obj, Pressure, dp)
            viscliquid = interp1(obj.data2.P(1,:),obj.data2.mu(1,:), Pressure*10^5);
            viscgas = interp1(obj.data2.P(2,:), obj.data2.mu(2,:), Pressure*10^5-dp );
            checkNaN(viscliquid,'Liquid viscosity', [], Pressure);
            checkNaN(viscgas,'Gas viscosity', [], Pressure-dp);
        end
        function [kliquid, kgas] = Bk(obj, Pressure, dp)
            kliquid = interp1(obj.data2.P(1,:),obj.data2.k(1,:), Pressure*10^5);
            kgas = interp1(obj.data2.P(2,:), obj.data2.k(2,:), Pressure*10^5-dp );
            checkNaN(kliquid,'Liquid k', [], Pressure);
            checkNaN(kgas,'Gas k', [], Pressure-dp);
        end
        function [densityliquid, densitygas] = Bdensity(obj, Pressure, dp)
            densityliquid = interp1(obj.data2.P(1,:),obj.data2.Rho(1,:), Pressure*10^5)*obj.data.Mw;
            densitygas = interp1(obj.data2.P(2,:), obj.data2.Rho(2,:), Pressure*10^5-dp )*obj.data.Mw;
            checkNaN(densityliquid,'Liquid density', [], Pressure);
            checkNaN(densitygas,'Gas density', [], Pressure-dp);
        end
        function [cpliquid, cpgas] = Bcp(obj, Pressure, dp)
            cpliquid = interp1(obj.data2.P(1,:),obj.data2.Cp(1,:), Pressure*10^5)/obj.data.Mw;
            cpgas = interp1(obj.data2.P(2,:), obj.data2.Cp(2,:), Pressure*10^5-dp)/obj.data.Mw;
            checkNaN(cpliquid,'Liquid Cp', [], Pressure);
            checkNaN(cpgas,'Gas Cp', [], Pressure-dp);
        end
        
        %function T_boiling = boiling(~, Pressure)
        %    p_table = (1:0.25:12);
        %    t_table = [20.324,21.102,21.775,22.371,22.91,23.404,23.86,24.285,24.683,25.059,25.415,25.753,26.076,26.386,26.682,26.968,27.243,27.508,27.765,28.014,28.255,28.489,28.717,28.939,29.154,29.365,29.57,29.771,29.967,30.158,30.346,30.529,30.709,30.885,31.057,31.227,31.393,31.555,31.715,31.872,32.027,32.178,32.327,32.473,32.616];
        %    T_boiling = interp1(p_table, t_table, Pressure,"linear");
        %end
        %function hvap = hvap(~, Temperature)
        %    T_table = [20.324,21.775,22.91,23.86,24.683,25.415,26.076,26.682,27.243,27.765,28.255,28.717,29.154,29.57,29.967,30.346,30.709,31.057,31.393,31.715,32.027,32.327,32.616];
        %    h_table = [448.917,440.759,431.955,422.756,413.266,403.529,393.534,383.292,372.773,361.957,350.81,339.31,327.4,314.99,302.04,288.43,274.04,258.69,242.13,224,203.7,180.2,151.34]*1000;
        %    if(Temperature<373.946+273.15)
        %        hvap = interp1(T_table,h_table,Temperature,"linear");
        %    else
        %        hvap =0;
        %    end
        %end
    end 
end