i = 0.1;
j = 1;
diam=0.1;
radiative_loss = zeros(5,50);
hr = 0.1:0.1:5;
while i<=5
    cone1 = cone(diam,diam*i,1*diam,1000,0.85);
    radiative_loss(1,j) = cone1.Radiation_Loss/(pi*cone1.Diameter/2*sqrt(cone1.Diameter^2/4+cone1.Length^2)+pi*(cone1.Diameter^2-cone1.ApertureDiameter^2)/4);
    cone2 = cone(diam, diam*i, 0.8*diam, 1000, 0.85);
    radiative_loss(2,j) = cone2.Radiation_Loss/(pi*cone2.Diameter/2*sqrt(cone2.Diameter^2/4+cone2.Length^2)+pi*(cone2.Diameter^2-cone2.ApertureDiameter^2)/4);
    cone3 = cone(diam,diam*i,0.6*diam,1000,0.85);
    radiative_loss(3,j) = cone3.Radiation_Loss/(pi*cone3.Diameter/2*sqrt(cone3.Diameter^2/4+cone3.Length^2)+pi*(cone3.Diameter^2-cone3.ApertureDiameter^2)/4);
    cone4 = cone(diam, diam*i, 0.4*diam, 1000, 0.85);
    radiative_loss(4,j) = cone4.Radiation_Loss/(pi*cone4.Diameter/2*sqrt(cone4.Diameter^2/4+cone4.Length^2)+pi*(cone4.Diameter^2-cone4.ApertureDiameter^2)/4);
    cone5 = cone(diam, diam*i, 0.2*diam, 1000, 0.85);
    radiative_loss(5,j) = cone5.Radiation_Loss/(pi*cone5.Diameter/2*sqrt(cone5.Diameter^2/4+cone5.Length^2)+pi*(cone5.Diameter^2-cone5.ApertureDiameter^2)/4);
 
    i = i+0.1;
    j= j+1;
end

    radiative_loss = radiative_loss/radiative_loss(1,1);
hold on
grid on 
title('Normalized Radiative Heat Loss')
xlabel('Cone Height/Diameter Ratio')
ylabel('Radiative Loss (W/m^2)')
plot(hr, radiative_loss(1,:),hr, radiative_loss(2,:),hr , radiative_loss(3,:), hr, radiative_loss(4,:), hr, radiative_loss(5,:));
legend({'AR = 1','AR = 0.8', 'AR = 0.6', 'AR = 0.4', 'AR = 0.2'},'Location','northeast')
hold off