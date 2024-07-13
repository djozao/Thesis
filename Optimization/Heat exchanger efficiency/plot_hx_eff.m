% Define the parameter ranges for h and w
h_range = 0.001:0.0005:0.003;  % Example range for h
w_range = 0.0005:0.0005:0.005;  % Example range for w
n_pipes = 1;
% Create a meshgrid for h and w
[H, W] = meshgrid(h_range, w_range);

% Initialize the efficiency matrix
efficiency = zeros(size(H));

% Compute the efficiency for each (h, w) pair
for i = 1:size(H, 1)
    for j = 1:size(W, 2)
        h = H(i, j);
        w = W(i, j);
        efficiency(i, j) = linear_duct_efficiency(5, h, w, 0.01, 251e-6, n_pipes, fluid, 700);
    end
end

% Create a 3D plot
figure;
surf(H, W, efficiency*100);

% Add labels
xlabel('Height (h) [mm]');
ylabel('Width (w) [mm]');
zlabel('Heat exchanger efficiency (%)');

% Add a color bar
colorbar;
