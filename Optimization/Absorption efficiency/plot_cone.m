function plot_cone(L, D, A, n_1, n_3, Matrix)

gamma = atan(D/(2*L));
L_D = L/cos(gamma);

x_apert = [0,0];
y_apert = [A/2, -A/2];

x_lat = zeros(1, n_3+1);
y_lat = zeros(1, n_3+1);
for i = 1:n_3+1
    x_lat(1,i) = L-L_D/n_3*(i-1)*cos(gamma);
    y_lat(1,i) = L_D/n_3*(i-1)*sin(gamma);
end

x_lat_sym = x_lat;
y_lat_sym = -y_lat;

x_base = zeros(1,n_1+1);
y_base = zeros(1,n_1+1);
for i = 1:n_1+1
    x_base(1,i) = 0;
    y_base(1,i) = D/2-(D-A)/2/n_1*(i-1);
end

% Normalize values to map to the colormap
min_value = min(Matrix);
max_value = max(Matrix);

% Create a figure to plot the data
figure;

cmap = colormap('jet'); % You can use a different colormap if needed

% Plot the first line segment with color based on the colormap
%color_segment = cmap(round((Matrix(1)-min_value)/(max_value-min_value) * (size(cmap, 1) - 1) + 1), :);
%plot(x_apert, y_apert,'Color', color_segment1, 'LineWidth', 2);

hold on
for i = 1:n_3
    color_segment = cmap(round((Matrix(n_1+i)-min_value)/(max_value-min_value) * (size(cmap, 1) - 1) + 1), :);
    plot([x_lat(i), x_lat(i+1)], [y_lat(i), y_lat(i+1)], 'Color', color_segment, 'LineWidth', 2 )
    plot([x_lat(i), x_lat(i+1)], [-y_lat(i), -y_lat(i+1)], 'Color', color_segment, 'LineWidth', 2 )
end
for i = 1:n_1
    color_segment = cmap(round((Matrix(i)-min_value)/(max_value-min_value) * (size(cmap, 1) - 1) + 1), :);
    plot([x_base(i), x_base(i+1)], [y_base(i), y_base(i+1)], 'Color', color_segment, 'LineWidth', 2 )
    plot([x_base(i), x_base(i+1)], [-y_base(i), -y_base(i+1)], 'Color', color_segment, 'LineWidth', 2 )
end

% Add a colorbar to indicate the value scale
c = colorbar;
c.Label.String = 'Temperature (K)';
caxis([min_value,max_value]);

% Add tick marks at minimum and maximum values
%ticks = [100*min_value, c.Ticks, 100*max_value];
%c.Ticks = ticks;


hold off

end