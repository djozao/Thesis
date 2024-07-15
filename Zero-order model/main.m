% Setup Run time
step = 60;  % step in seconds
time = 60*360; % run time in seconds
time_table = (step:step:time);  % time data
i = 1; % accounts for the step increase



% Setup RAC Geometry and properties
CR = 0; % cone ratio: 0 - cylinder 1- cone
diameter = 0.04;%m
aperture = 0.04;%m
length = 0.1; %m
abs = 0.1; % absorptivity
emi = 0.1; % emissivity
mass = 1; %kg
c_p = 900; % J/(K*kg)
dif = 1; % 0 means specular and 1 for diffuse reflection
cavity = geometry(CR, length, diameter, aperture, emi);

% Beam properties
P_laser = 100; % W
alpha = 5; % half-beam angle in degrees
o = aperture/(2*tan(alpha)); % origin of beam in meters


% ducts properties
numberpipes = 8;
L = 0.1; % duct length in meters
A = 0.5*1*10^-6;  % area of each duct in m^2
P=0.5/1000*2+1/1000*2; % perimeter of each duct in m^2

% Propellant
T_pi = 293.15;  % inlet temperature in K
T_p = zeros(time/step,1);  % saves temperature in time
dp = zeros(time/step,1);   %
P_i = 5;   % initial pressure in bar 
propellant = nitrogen(200:50:1000, 3:0.25:5);   % define propellant and acquire properties from NIST
mflow  = 250e-6/numberpipes;  % mass flow in each duct in kg/s

% Variables ambient and setting 
T_amb = 293.15;  % ambient temperature in K
T_w = T_amb;  % initial RAC temperature in K

Tw = zeros(time/step,1);  % saves RAC wall temperature transient data
Radloss = zeros(time/step,1); % saves heat loss by inner radiation transient data
Convloss = zeros(time/step,1); % saves heat loss by inner convection transient data
Radloss_outer = zeros(time/step,1); % saves heat loss by outer radiation transient data
Convloss_outer = zeros(time/step,1); % saves heat loss by outer convection transient data
Q_propel = zeros(time/step,1); % saves transient data of the heat transferred to the propelant

%Tins = zeros(time/step,1);
%Condloss = zeros(time/step,1);
%T_ins = T_amb;

P_abs = power_absorbed_function(P_laser, alpha, length, diameter, aperture, CR, abs, dif);




%%  Set up the transient conditions; 
while i<=time/step
    
    % Calculate the final propellant temperature, heat transfered to each duct, pressure
    % drop and the sections where the propellant is single phase and two
    % phase flow
    % For spiral ducts use the spiral_duct function and add thee diameter
    [T_pf, Q_prop, dp1, sections] = linear_duct(T_pi, T_w, P_i, A, P, L, mflow, propellant);
    % Calculate the inner radiation and convection loss
    Radloss(i) = cavity.radiationloss(T_w,T_amb);
    Convloss(i) = 0; %cavity.convectionloss(T_w,T_amb);
    % Calculate the overall propellant heat flux
    Q_propel(i) = Q_prop*numberpipes;
    % Register the final propellant temperature, pressure drop and RAC
    % temperature
    T_p(i)= T_pf;
    dp(i) = dp1*numberpipes;
    Tw(i) = T_w;
    
    % Insulation calculation for the used in Leenders with its respective
    % thermal conductivity formula from Takken's thesis
    %Tins(i) = T_ins;
    %[T_ins, Radloss_outer(i), Convloss_outer(i), Condloss(i)] = T_insulation_calculation(cavity.Length,cavity.Diameter,cavity.ApertureDiameter,thickness_ins,T_w,T_amb,0.1);
    

    % Calculate outer losses (radiation and convection) if CR != 0 || CR
    % !=1 the convection loss could be estimated with an averaged approach
    % for example if CR = 0.1 use 0.1*convloss_cone+0.9*convloss_cyl
    Radloss_outer(i) = 0; %(T_w^4-T_amb^4)*outer_radlossfrac_cylinder(length, diameter, aperture, emi);
    Convloss_outer(i) = 0;% outer_convloss_cylinder(cavity.Length,cavity.Diameter,cavity.ApertureDiameter,T_w,T_amb);
    
    % Update the RAC temperature considering the heat flows calculated
    T_w = T_w+step*(P_abs-Q_propel(i)-Radloss(i)-Convloss(i)-Radloss_outer(i)-Convloss_outer(i))/(mass*c_p)

    i = i+1
end











