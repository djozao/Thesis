classdef cone
    properties (SetAccess = public)
        Diameter (1,1) double
        Length (1,1) double
        ApertureDiameter (1,1) double
        Emissivity (1,1) double
    end
    methods
        function obj = cone(diameter, apertureDiameter, length, emi)
            if nargin > 0
                obj.Diameter = diameter;
                obj.Length = length;
                obj.ApertureDiameter = apertureDiameter;
                obj.Emissivity = emi;
            end
        end
        function rad_loss = radiationloss(obj, Temperature, T_amb)
               rad_loss = (Temperature^4-T_amb^4)*radiationlossfrac_cone(obj.Length, obj.Diameter, obj.ApertureDiameter, obj.Emissivity);
        end 
        function conv_loss = convectionloss(obj, Temperature, T_amb)
            cyl = cylinder (obj.Diameter,obj.Length, obj.ApertureDiameter, Temperature, obj.Emissivity);
            conv_loss = cyl.Convection_Loss(Temperature,T_amb)*0.75;
        end
    end
end