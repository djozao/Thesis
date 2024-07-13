classdef nitrogen
    properties (SetAccess = public)
        data
        data2
    end
    methods
        function obj = nitrogen(Temperature, Pressure)
            if nargin > 0
                 [obj.data,obj.data2] = nistdata('N2',Temperature,Pressure);
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
    end 
end