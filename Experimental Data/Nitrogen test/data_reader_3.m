data = tdmsread("data_3.tdms")

%% mass flow
time_mf = data{1,1}.time;
massflow = data{1,1}.('Feed system mass flow {Mass flow} [mgps]');
plot(time_mf, massflow)

% Add labels and title
xlabel('Time (s)');
ylabel('Mass Flow Rate (mg/s)');

% Display grid
grid on;


%%
% Define time intervals
intervals = [51, 61; 110, 120; 170, 180; 231, 241; 290, 300]; % Define your intervals here

% Initialize arrays to store averaged mass flow
averaged_mass_flow = zeros(size(intervals, 1), 1);
std_mass_flow = zeros(size(intervals,1), 1);
% Loop through each time interval
for i = 1:size(intervals, 1)
    % Find indices corresponding to the current interval
    indices = find(time_mf >= intervals(i, 1) & time_mf <= intervals(i, 2));
    
    % Calculate average mass flow within the interval
    averaged_mass_flow(i) = mean(massflow(indices));
    std_mass_flow(i) = std(massflow(indices));
end

% Display averaged mass flow for each interval
for i = 1:size(intervals, 1)
    fprintf('Average mass flow between %d and %d: %.2f\n', intervals(i, 1), intervals(i, 2), averaged_mass_flow(i));
end

%% Pressure

filename = '20240325 THIRD\pressure_3.txt'; % Replace 'your_file.txt' with the actual filename
fid = fopen(filename, 'r');
data_p = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

% Extract number of rows
num_rows = numel(data_p{1});

% Initialize arrays to store extracted values
time_pressure = zeros(num_rows, 1);
pressure1 = zeros(num_rows, 1);
pressure2 = zeros(num_rows, 1);

% Loop through each line of data
for i = 1:num_rows
    % Split the line into individual values
    values = strsplit(data_p{1}{i}, '\t');
    
    % Extract time, pressure1, and pressure2 values
    % Time
    time_str = values{1}; % Time string
    time_pressure(i) = str2double(time_str(1:end-4)); % Extract time without decimals
    
    % Pressure1
    pressure1(i) = 1+str2double(values{2})/1000; % Convert string to numeric
    
    % Pressure2
    pressure2(i) = 1+str2double(values{3})/1000; % Convert string to numeric
end

% Subtract the first entry of time and divide by 1000 to convert to seconds
time_pressure = (time_pressure - time_pressure(1)) / 1000;

% Now time will start from 0 and will be in seconds

% Define smoothing window size
window_size = 10; % Adjust as needed

% Smooth the pressure2 data
smoothed_pressure2 = smooth(pressure2, window_size);
figure
% Plot pressure2 vs. time
plot(time_pressure, pressure2, 'b-', 'DisplayName', 'Original Pressure2'); % Plot original pressure2 in blue
hold on;

% Plot smoothed pressure2 vs. time
%plot(time_pressure, smoothed_pressure2, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed Pressure2'); % Plot smoothed pressure2 in red with thicker line

% Add labels and title
xlabel('Time (s)');
ylabel('Absolut Inlet Pressure (bar)');

% Display legend
%legend();

% Display grid
grid on;

% Show the plot


% Define time intervals (in seconds)
intervals = [35, 61; 95, 121; 154, 181; 215, 241 ; 275, 302]; % Adjust as needed

% Initialize arrays to store average and standard deviation values
avg_pressure = zeros(size(intervals, 1), 1);
std_pressure = zeros(size(intervals, 1), 1);

% Loop through each time interval
for i = 1:size(intervals, 1)
    % Find indices corresponding to the current interval
    interval_indices = find(time_pressure >= intervals(i, 1) & time_pressure <= intervals(i, 2));
    
    % Extract pressure values within the interval
    pressure_interval = pressure2(interval_indices);
    
    % Calculate average and standard deviation
    avg_pressure(i) = mean(pressure_interval);
    std_pressure(i) = std(pressure_interval);
end

% Display results
for i = 1:size(intervals, 1)
    fprintf('Interval %d: Average Pressure = %.2f, Standard Deviation = %.2f\n', i, avg_pressure(i), std_pressure(i));
end


%%

time_load = data{1,2}.time;
load = data{1,2}.('V_LC [V]');
plot(time_load, load)

% Smooth the data using a moving average filter
windowSize = 100; % Adjust this value to change the smoothing effect
smoothedY = movmean(load, windowSize);
figure
% Plot original and smoothed data
plot(time_load, load, 'b', 'LineWidth', 1.5);
hold on;
plot(time_load, smoothedY, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Load Cell Output Voltage (V)');
legend('Original', 'Smoothed');
grid on;
hold off;

%%

% Define time intervals (x and y values)
time_intervals_no_load = [0, 30; 72, 90; 132, 150; 191, 209; 253, 270];
time_intervals_load = [36, 61; 94, 121 ; 155, 181 ; 215, 242 ; 276, 302];

% Initialize a cell array to store segmented load data
segmented_load_data_no_load = cell(size(time_intervals_no_load, 1), 1);
segmented_load_data_load = cell(size(time_intervals_load, 1), 1);

averaged_load_data_no_load = zeros(size(time_intervals_no_load,1),1);
averaged_load_data_load = zeros(size(time_intervals_load,1),1);

std_load_data_no_load = zeros(size(time_intervals_load,1),1);
std_load_data_load = zeros(size(time_intervals_load,1),1);


% Iterate through each time interval
for i = 1:size(time_intervals_no_load, 1)
    % Extract time and load data within the current time interval
    x = time_intervals_no_load(i, 1);
    y = time_intervals_no_load(i, 2);
    segment_indices = time_load >= x & time_load <= y;
    segmented_load_data_no_load{i} = load(segment_indices);
    averaged_load_data_no_load(i) = mean(segmented_load_data_no_load{i});
    std_load_data_no_load(i) = std(segmented_load_data_no_load{i});
end

% Iterate through each time interval
for i = 1:size(time_intervals_load, 1)
    % Extract time and load data within the current time interval
    x = time_intervals_load(i, 1);
    y = time_intervals_load(i, 2);
    segment_indices = time_load >= x & time_load <= y;
    segmented_load_data_load{i} = load(segment_indices);
    averaged_load_data_load(i) = mean(segmented_load_data_load{i});
    std_load_data_load(i) = std(segmented_load_data_load{i});
end

differences = averaged_load_data_load-averaged_load_data_no_load;

for i = 1:size(time_intervals_no_load, 1)
    fprintf('Interval %d: Average non-Load = %.3f, Standard Deviation = %.3f\n', i, averaged_load_data_no_load(i), std_load_data_no_load(i));
end
for i = 1:size(time_intervals_load, 1)
    fprintf('Interval %d: Average Load = %.3f, Standard Deviation = %.3f\n', i, averaged_load_data_load(i), std_load_data_load(i));
end

for i = 1:size(differences)
    fprintf('Interval %d: Average difference = %.3f \n', i, differences(i));
end

%%

time_load = data{1,2}.time;
load_data = data{1,2}.('V_LC [V]');
plot(time_load, load_data, 'b', 'LineWidth', 1.5);
hold on;

% Smooth the data using a moving average filter
windowSize = 100; % Adjust this value to change the smoothing effect
smoothedY = movmean(load_data, windowSize);
plot(time_load, smoothedY, 'r', 'LineWidth', 2);

% Define time intervals
time_intervals_no_load = [0, 30; 72, 90; 132, 150; 191, 209; 253, 270];
time_intervals_load = [36, 61; 94, 121 ; 155, 181 ; 215, 242 ; 276, 302];

% Plot vertical lines for time_intervals_no_load
for i = 1:size(time_intervals_no_load, 1)
    line([time_intervals_no_load(i, 1), time_intervals_no_load(i, 1)], ylim, 'Color', 'g', 'LineStyle', '--');
    line([time_intervals_no_load(i, 2), time_intervals_no_load(i, 2)], ylim, 'Color', 'g', 'LineStyle', '--');
end

% Plot vertical lines for time_intervals_load
for i = 1:size(time_intervals_load, 1)
    line([time_intervals_load(i, 1), time_intervals_load(i, 1)], ylim, 'Color', 'm', 'LineStyle', '--');
    line([time_intervals_load(i, 2), time_intervals_load(i, 2)], ylim, 'Color', 'm', 'LineStyle', '--');
end

xlabel('Time');
ylabel('Load');
title('Original vs. Smoothed Data with Time Intervals');
legend('Original', 'Smoothed', 'No Load Intervals', 'Load Intervals');
grid on;
hold off;

%% Calculation


mN_A_stiff = 0.7926;
mN_A_free = 0.8253;
d_thruster = 187.5;
d_loadcell = 175;
V_mN_loadcell = 0.043353720693170;

ratio_free_stiff = mN_A_free/mN_A_stiff;
ratio_Treal_Tmeasured = d_loadcell/d_thruster;

measured_thrusts = differences/V_mN_loadcell
Thrust_real_stiff = measured_thrusts*ratio_Treal_Tmeasured
Thrust_real_free = Thrust_real_stiff*ratio_free_stiff 



%% Expected Thrust
gamma_nitrogen = 1.40;
Vanden =  sqrt(gamma_nitrogen)*(2/(gamma_nitrogen+1))^((1+gamma_nitrogen)/(2*(gamma_nitrogen-1)));
nozzle_diameter = 0.5*10^-3;
A_t = pi*(nozzle_diameter)^2/4;
A_e = A_t;
T_c = 273.15+14.918;
R_nitrogen = 296.80;
p_a = 1.066;

dp = 0.29;
p_c = avg_pressure-dp;

p_t = p_c*(2/(gamma_nitrogen+1))^(gamma_nitrogen/(gamma_nitrogen-1));
p_e = p_t;

mf_expected = Vanden*A_t*p_c*10^5/sqrt(T_c*R_nitrogen) %kg/s

mf_adjusted = (averaged_mass_flow-(0.1257.*(avg_pressure-p_a).^2+3.9322*(avg_pressure-p_a)-0.006  ) )*10^-6 %kg/s

U_e = sqrt(2*(gamma_nitrogen)/(gamma_nitrogen-1)*R_nitrogen*T_c*(1-(p_e./p_c).^((gamma_nitrogen-1)/gamma_nitrogen)))
U_eq = U_e + (p_e-p_a)*10^5./mf_expected*A_e % m/s

U_eq_measured_mf = U_e + (p_e-p_a)*10^5./(averaged_mass_flow*10^-6)*A_e % m/s

U_eq_adjusted_mf = U_e + (p_e-p_a)*10^5./(mf_adjusted)*A_e


F_theory = mf_expected.*U_eq


F_theory_measured_mf = (averaged_mass_flow*10^-6).*U_eq_measured_mf *10^3


F_theory_adjusted_mf = mf_adjusted.*U_eq_adjusted_mf*10^3

T_t = mean(U_e)^2/(gamma_nitrogen*R_nitrogen);

coef = Thrust_real_free./F_theory_adjusted_mf

C_f_real = mean(Thrust_real_free/10^3./p_c/10^5/A_t/(mean(mf_adjusted./mf_expected)));
C_f_theory = Vanden*sqrt(2*gamma_nitrogen/(gamma_nitrogen-1)*(1-(mean(p_e./p_c))^((gamma_nitrogen-1)/gamma_nitrogen))) +mean((p_e-p_a)./p_c);

nozzle_flow_quality = C_f_real/C_f_theory

c_star_real = mean(p_c*10^5*A_t./mf_expected)
c_star_theory = 1/Vanden*sqrt(R_nitrogen*T_c)

%F_maybe = 1/2.*mf_adjusted.*U_e + (p_e-p_a)*10^5*A_e;

%coef_maybe = Thrust_real_free./F_maybe/10^3;

%% Reynolds number

viscosity_throat = 1.5022*10^-5;

Re_throat = 4*mean(mf_adjusted)/viscosity_throat/(pi*nozzle_diameter)


%C_f_vis= 17.6*exp((0.0032*A_e/A_t))./(sqrt(0.773.*Re_throat))

%%
mf_mean = mean(mf_adjusted)
viscosity = 1.7381*10^-5;
density = 4.6823;

Re_first = 4*mf_mean/(viscosity*pi*3.3*10^-3)
Re_second = 4*mf_mean/8/(viscosity*3*10^-3)
Re_third = 4*mf_mean/(viscosity*pi*1*10^-3)

f_first = 1/((1.82*log10(Re_third)-1.64)^2)
f_second = 64/Re_second
f_third = 1/((1.82*log10(Re_third)-1.64)^2)

loss_p_first = f_first * 0.07/(0.0033)*1/2*(mf_mean/(pi*(3.3*10^-3)^2/4))^2/density
loss_p_second = f_second * 0.08/(4*1*0.5*10^-6/(3*10^-3))*1/2*(mf_mean/8/(1*0.5*10^-6))^2/density*8
loss_p_third = f_third*0.25/(0.001)*1/2*(mf_mean/(pi*(0.001)^2/4))^2/density

loss = (loss_p_third+loss_p_second+loss_p_first)/10^5
%%
mean(mf_expected*10^6)
mean(mf_adjusted*10^6)
mean(mf_adjusted*10^6)/mean(mf_expected*10^6)
Thrust_real_free./(F_theory)./(mf_adjusted./mf_expected)*10^-3
nozzle_flow_quality
mean((F_theory))


