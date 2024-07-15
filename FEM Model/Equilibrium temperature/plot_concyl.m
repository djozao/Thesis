function plot_concyl(L, D, A, CR, n1, n2, n3, Matrix)

gamma = atan(D/(2*L*CR));
L_D = L*CR/cos(gamma);


x_lat = zeros(1, n2+1);
y_lat = zeros(1, n2+1);
for i = 1:n2+1
    x_lat(1,i) = L*(1-CR)-L*(1-CR)/n2*(i-1);
    y_lat(1,i) = D/2;
end

x_base = zeros(1,n1+1);
y_base = zeros(1,n1+1);
for i = 1:n1+1
    x_base(1,i) = 0;
    y_base(1,i) = D/2-(D-A)/2/n1*(i-1);
end

x_bottom = zeros(1, n3+1);
y_bottom = zeros(1, n3+1);
for i = 1:n3+1
    x_bottom(1,i) = L-L_D/n3*(i-1)*cos(gamma);;
    y_bottom(1,i) = D/2/n3*(i-1);
end

% Normalize values to map to the colormap
min_value = min(Matrix);
max_value = max(Matrix);

% Create a figure to plot the data
figure;

cmap = colormap('jet'); % You can use a different colormap if needed

% Plot the first line segment with color based on the colormap
%color_segment = cmap(round((Heat(1)-min_value)/(max_value-min_value) * (size(cmap, 1) - 1) + 1), :);
%plot(x_apert, y_apert,'Color', color_segment1, 'LineWidth', 2);

hold on
for i = 1:n1
    color_segment = cmap(round((Matrix(i)-min_value)/(max_value-min_value) * (size(cmap, 1) - 1) + 1), :);
    plot([x_base(i), x_base(i+1)], [y_base(i), y_base(i+1)], 'Color', color_segment, 'LineWidth', 2 )
    plot([x_base(i), x_base(i+1)], [-y_base(i), -y_base(i+1)], 'Color', color_segment, 'LineWidth', 2 )
end
for i = 1:n2
    color_segment = cmap(round((Matrix(n1+i)-min_value)/(max_value-min_value) * (size(cmap, 1) - 1) + 1), :);
    plot([x_lat(i), x_lat(i+1)], [y_lat(i), y_lat(i+1)], 'Color', color_segment, 'LineWidth', 2 )
    plot([x_lat(i), x_lat(i+1)], [-y_lat(i), -y_lat(i+1)], 'Color', color_segment, 'LineWidth', 2 )
end
for i = 1:n3
    color_segment = cmap(round((Matrix(n1+n2+i)-min_value)/(max_value-min_value) * (size(cmap, 1) - 1) + 1), :);
    plot([x_bottom(i), x_bottom(i+1)], [y_bottom(i), y_bottom(i+1)], 'Color', color_segment, 'LineWidth', 2 )
    plot([x_bottom(i), x_bottom(i+1)], [-y_bottom(i), -y_bottom(i+1)], 'Color', color_segment, 'LineWidth', 2 )
end

% Add a colorbar to indicate the value scale
c = colorbar;
c.Label.String = 'Temperature (K)';
caxis([min_value,max_value]);

% Add tick marks at minimum and maximum values
%ticks = [100*min_value, c.Ticks, 100*max_value];
%c.Ticks = ticks;


hold off
%%


end