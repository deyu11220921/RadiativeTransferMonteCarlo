%% Parameter Sweep and Calibration Script
clear; clc; close all;

%% 1. Basic Settings
% Complete geometry settings to prevent initialization errors
geometry.dimension = 3; 
geometry.frame = 'spherical'; 
geometry.bnd(1) = struct('dir',3,'val',0);

freq = 1e6;        % Frequency: 1 MHz
V_paper = 450;    % S-wave velocity from reference paper (Granite)
acoustics = true; 
ac_func = 'exp';   % Autocorrelation function: exponential

% Three target mean free paths (ls) based on paper (L/ls = 1, 3, 10)
% Assuming sample length L = 10mm
target_ls_list = [10e-3, 3.3e-3, 1e-3]; 
target_names = {'Weak Scattering (L/ls=1)', 'Medium (L/ls=3)', 'Strong (L/ls=10)'};

%% 2. Define Scanning Grid
N = 25; 
variance_values = linspace(0.05, 0.2, N);    % Range for Epsilon (Variance)
correlation_values = linspace(5e-4, 2.0e-3, N); % Range for Correlation Length (a)
Results_ls = zeros(N, N);

%% 3. Execution Loop (Grid Search)
fprintf('Calculating % d x %d parameter sets...\n', N, N);

for i = 1:N 
    for j = 1:N 
        current_var = variance_values(i);
        current_lc = correlation_values(j);
        
        % Initialize MaterialClass with current sweep parameters
        mat = MaterialClass(geometry, freq, acoustics, V_paper, ...
                            [current_var current_var], ...
                            -0.5, ... 
                            ac_func, current_lc);
        
        % Calculate the scattering covariance matrix/tensor
        mat = prepareSigma(mat, geometry.dimension); 
        
        % Store the resulting Mean Free Path (ls)
        Results_ls(i,j) = mat.meanFreePath; 
    end
end

%% 4. Automatic Best-Fit Parameter Identification
fprintf('\n================ Calibration Results ================\n');
fprintf('Copy these parameters into mainAcoustics for simulation:\n\n');

for t = 1:length(target_ls_list)
    t_val = target_ls_list(t);
    
    % Find the value in the results matrix closest to the target ls
    diff_matrix = abs(Results_ls - t_val);
    [min_val, idx] = min(diff_matrix(:));
    [row, col] = ind2sub(size(diff_matrix), idx);
    
    best_epsilon = variance_values(row);
    best_a = correlation_values(col);
    calc_ls = Results_ls(row, col);
    
    fprintf('Target Case: %s\n', target_names{t});
    fprintf('  -> Target ls = %.4f mm\n', t_val*1000);
    fprintf('  -> Best Match: Epsilon (Var) = %.4f, a (Corr Length) = %.4f mm\n', best_epsilon, best_a*1000);
    fprintf('  -> Calculated ls = %.4f mm (Error %.2f%%)\n', calc_ls*1000, abs(calc_ls-t_val)/t_val*100);
    fprintf('------------------------------------------\n');
end

%% 5. Visualization (3D Trend Analysis & Target Calibration)
% Based on teacher's requirements: observe ls variation in 3D 
[X_grid, Y_grid] = meshgrid(correlation_values, variance_values);

figure('Color', 'w', 'Name', 'Mean Free Path Calibration');
hold on; grid on;

% 1. Create the 3D Surface Plot as requested by the teacher 
% X: Correlation Length (mm), Y: Epsilon (Variance), Z: Mean Free Path (mm)
s = surf(X_grid*1000, Y_grid, Results_ls*1000, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
colormap(jet);
colorbar;
xlabel('Correlation Length a (mm)');
ylabel('Epsilon (Variance)');
zlabel('Mean Free Path l_s (mm)');
title('3D Mean Free Path Variation & Target Selection');

% 2. Draw 3D Contour lines for target scattering regimes [cite: 355-365]
% These lines represent the "Solution Lines" for L/ls = 1, 3, 10
L_ls_ratios = [1, 3, 10];
[C, h] = contour3(X_grid*1000, Y_grid, Results_ls*1000, target_ls_list*1000, ...
                  'r-', 'LineWidth', 3);
clabel(C, h, 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');

% 3. Mark the calibrated "Best Match" points on the surface
marker_colors = {'g', 'm', 'c'}; 
for t = 1:length(target_ls_list)
    t_val = target_ls_list(t);
    diff_matrix = abs(Results_ls - t_val);
    [~, idx] = min(diff_matrix(:));
    [row, col] = ind2sub(size(diff_matrix), idx);
    
    % Plot point exactly on the surface
    plot3(correlation_values(col)*1000, variance_values(row), Results_ls(row,col)*1000, ...
          'o', 'MarkerSize', 10, 'MarkerFaceColor', marker_colors{t}, ...
          'MarkerEdgeColor', 'k', 'LineWidth', 2);
    
    % Add label for the scattering regime
    text(correlation_values(col)*1000, variance_values(row), Results_ls(row,col)*1000 + 0.5, ...
         sprintf(' L/l_s = %d', L_ls_ratios(t)), 'FontSize', 11, 'FontWeight', 'bold');
end

% Set initial view to show the 3D trend 
view(-35, 30); 

% Optional: Use 'view(0, 90)' later to check parameter alignment from top-down