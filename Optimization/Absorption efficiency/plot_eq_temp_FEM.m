function plot_eq_temp_FEM()
    % Define optimal solution
    x_opt = [5, 0.05, 0.04, 0.03999, 0.15, 0.05, 0.004]; % Example optimal solution
    P = 20;
    CR = 0;
    n_1 = 10;
    n_2 = 10;
    n_3 = 10;
    k = 237;
    % Define intervals for L and D
    L_interval = 0.05:0.01:0.1;
    D_interval = 0.04:0.01:0.1;
    %Angle_interval = 5:1:20;

    % Generate grid of L and D values
    [L_grid, D_grid] = meshgrid(L_interval, D_interval);
    %[Angle_grid, L_grid2] = meshgrid(Angle_interval, L_interval);

    % Calculate equilibrium temperature for each combination of L and D
    T_eq = zeros(size(L_grid));
    eff = zeros(size(L_grid));
    for i = 1:numel(L_grid)
        if (D_grid(i)>= x_opt(4))
        [T_eq(i), eff(i)] = CalcEqTemp_FEM_Geo(CR, P, x_opt(1), L_grid(i), D_grid(i), x_opt(4) , n_1, n_2, n_3, x_opt(5), x_opt(6), x_opt(7), k, 50);
        else
            T_eq(i) = 50;
        end
        i;
    end

    % Plot 3D surface plot
    figure;
    surf(L_grid, D_grid, reshape(eff, size(L_grid)));
    xlabel('Length (m)');
    ylabel('Diameter (m)');
    zlabel('Equilibrium Temperature (K)');
    colorbar;

    max_column = max(T_eq);
    max_T_eq = max(max_column)


    % T_eq_2 = zeros(size(Angle_grid));
    % for i = 1:numel(Angle_grid)
    % 
    %     T_eq_2(i) = CalcEqTemp_FEM_Geo(CR, P, Angle_grid(i), L_grid2(i), x_opt(3), x_opt(4), n_1, n_2, n_3, x_opt(5), x_opt(6), x_opt(7), k, 50);
    % end
    % 
    % % Plot 3D surface plot
    % figure;
    % surf(Angle_grid, L_grid2, reshape(T_eq_2, size(Angle_grid)));
    % xlabel('Angle (º)');
    % ylabel('Length (m)');
    % zlabel('Equilibrium Temperature (K)');
    % colorbar;

end