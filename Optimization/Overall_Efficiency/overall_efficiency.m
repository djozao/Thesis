%function eff = overall_efficiency()
eff_m2 = zeros(6,7);
P_abs_m2 = zeros(6,7);
T_w_m2 = zeros(6,7);
T_p_m2 = zeros(6,7);
zero_m2 = zeros(6,7);
Q_prop_m2 = zeros(6,7);
j_l = 1;
j_D = 1;
for L = 0.05:0.025:0.1
    for D = 0.04:0.03:0.1
% Geometry
Cyl = 0;
%D = 0.04;
A = 0.039999;
%L = 0.08;
Abs = 0.15;
mass = 0.2;
alpha = 5;


stefan = 5.670374*10^-8;
inputData = [alpha, Abs, A/D, L/D];
% Convert the numeric array to a table if needed
inputTable = array2table(inputData, 'VariableNames', {'Angle', 'Emi', 'Aperture', 'Length'});
emi = 0.05;
%Gebhart_matrix = Gebhart_3D_Cone(L,D,A,emi, 1, 1);
Gebhart_matrix = Gebhart_3D_Cylinder(L,D,A,emi, 1, 1, 1);
Area = [pi*(D^2-A^2)/4, pi*D*L, pi*D^2/4];
%Area = [pi*(D^2-A^2)/4, pi*D/2*sqrt(L^2+(D/2)^2)];
Rad_frac = stefan*emi*(Area*Gebhart_matrix(2:4,1));

% Variables ambient 
T_amb = 293.15;
T_w = T_amb;
P_laser = 20;
load("C:\Users\Diogo\Desktop\THESIS\Optimzation Tool\models\dif_CC00.mat");
P_abs = P_laser*dif_CC00.predictFcn(inputTable)/100
if P_abs >P_laser || P_abs<0
    P_abs = P_laser;
end
%ducts
numberpipes = 8;
h = 0.003;
w = 0.0005;
A = h*w;
P= (h+w)*2;


%transient
step = 60;
time = 60*1000;
time_table = (step:step:time);
i = 1;

% propellant
T_pi = 293.15;
T_p = zeros(time/step,1);
dp = zeros(time/step,1);
P_i = 6;
propellant = nitrogen(100:50:1000, 1:0.5:7);
mflow = 2500e-6/numberpipes;



Tw = zeros(time/step,1);
Radloss = zeros(time/step,1);
Q_propel = zeros(time/step,1);

i = 1;

while i< 1000
    if i< 10
    Q_prop = 0;
    T_pf = T_pi;
    dp1 = 0;
    else
    [T_pf, Q_prop, dp1] = singlephase_linear_duct(T_pi, T_w, P_i, A, P, L, mflow, propellant);
    end
    
    Radloss(i) = Rad_frac*(T_w^4-T_amb^4);
    Q_propel(i) = Q_prop*numberpipes;
    T_p(i)= T_pf;
    dp(i) = dp1*numberpipes;
    Tw(i) = T_w;
    T_w = T_w+step*(P_abs-Q_propel(i)-Radloss(i))/(mass*1000);
    T_w
    T_pf
    %T_ins
    i = i+1;
end

Q_prop;
Rad_frac*(T_w^4-T_amb^4);
Q_prop*8+Rad_frac*(T_w^4-T_amb^4)-P_abs
eff = Q_prop*numberpipes/P_laser;
T_w;

%Z(1,j) = P_abs
%Z(2,j) =T_w
%Z(3,j) =T_pf
%Z(4,j) = eff
%Z(5,j) = Q_prop*8+Rad_frac*(T_w^4-T_amb^4)-P_abs
%Z(6,j) = Q_prop*8
%Z(7,j) = Rad_frac*(T_w^4-T_amb^4)
%j = j+1;
    eff_m2(j_l, j_D) = eff; 
    P_abs_m2(j_l, j_D) = P_abs;
    T_w_m2(j_l, j_D) = T_w;
    T_p_m2(j_l, j_D) = T_pf;
    zero_m2(j_l, j_D) = Q_prop*8+Rad_frac*(T_w^4-T_amb^4)-P_abs;
    Q_prop_m2(j_l,j_D) = Q_prop*8;
    j_D = j_D +1;
    
    end
    j_D = 1;
    j_l = j_l+1;
end


