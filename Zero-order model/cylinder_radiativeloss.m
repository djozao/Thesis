i = 0.1;
j = 1;
radiative_loss = zeros(5,50);
hr = 0.1:0.1:5;
while i<=5
    cyl = cylinder(0.1,0.1*i,0.1,1000,0.85);
    radiative_loss(1,j) = cyl.Radiation_Loss;%/(pi*cyl.Diameter^2/4+pi*cyl.Diameter*cyl.Length+pi/4*(cyl.Diameter^2-cyl.ApertureDiameter^2));
    cyl2 = cylinder(0.1,0.1*i,0.08,1000,0.85);
    radiative_loss(2,j) = cyl2.Radiation_Loss;%/(pi*cyl2.Diameter^2/4+pi*cyl2.Diameter*cyl2.Length+pi/4*(cyl2.Diameter^2-cyl2.ApertureDiameter^2));
    cyl3 = cylinder(0.1,0.1*i,0.06,1000,0.85);
    radiative_loss(3,j) = cyl3.Radiation_Loss;%/(pi*cyl3.Diameter^2/4+pi*cyl3.Diameter*cyl3.Length+pi/4*(cyl3.Diameter^2-cyl3.ApertureDiameter^2));
    cyl4 = cylinder(0.1,0.1*i,0.04,1000,0.85);
    radiative_loss(4,j) = cyl4.Radiation_Loss;%/(pi*cyl4.Diameter^2/4+pi*cyl4.Diameter*cyl4.Length+pi/4*(cyl4.Diameter^2-cyl4.ApertureDiameter^2));
    cyl5 = cylinder(0.1,0.1*i,0.02,1000,0.85);
    radiative_loss(5,j) = cyl5.Radiation_Loss;%/(pi*cyl5.Diameter^2/4+pi*cyl5.Diameter*cyl5.Length+pi/4*(cyl5.Diameter^2-cyl5.ApertureDiameter^2));
    i = i+0.1;
    j= j+1;
end

    %radiative_loss = radiative_loss/radiative_loss(1,1);
hold on
grid on 
title('Normalized Radiative Heat Loss')
xlabel('Cylinder Height/Diameter Ratio')
ylabel('Radiative Loss % (1/m^2)')
plot(hr, radiative_loss(1,:), hr, radiative_loss(2,:),hr, radiative_loss(3,:),hr, radiative_loss(4,:),hr, radiative_loss(5,:));
legend({'AR = 1','AR = 0.8', 'AR = 0.6', 'AR = 0.4', 'AR = 0.2'},'Location','northeast')
hold off
