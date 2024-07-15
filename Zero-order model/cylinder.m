classdef cylinder
    properties (SetAccess = public)
        Diameter (1,1) double
        Length (1,1) double
        ApertureDiameter (1,1) double
        Emissivity (1,1) double
    end
    methods
        function obj = cylinder(diameter, apertureDiameter, length, emi)
            if nargin > 0
                obj.Diameter = diameter;
                obj.Length = length;
                obj.ApertureDiameter = apertureDiameter;
                obj.Emissivity = emi;
            end
        end 
        function rad_loss = radiationloss(obj, Temperature, T_amb)
               rad_loss = (Temperature^4-T_amb^4)*radiationlossfrac_cylinder(obj.Length, obj.Diameter, obj.ApertureDiameter, obj.Emissivity);
        end 
        function conv_loss = convectionloss(obj, Temperature, T_amb)
                conv_loss = convectionloss_cylinder(obj.Length,obj.Diameter,obj.ApertureDiameter,Temperature, T_amb);
        end
        function outer_radloss = outer_radiationloss(obj,Temperature,T_amb,emi)
            outer_radloss = outer_radlossfrac_cylinder(obj.Length, obj.Diameter, obj.ApertureDiameter, emi)*(Temperature^4-T_amb^4); 
        end
    end
end



