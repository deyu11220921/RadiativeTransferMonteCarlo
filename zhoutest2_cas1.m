% mainAcoustic_ZHOU2021.m
% Reproduction of Figure 1 from Zhou et al. (2021)
% Case: Infinite Medium (Comparison of L = 1*ls, 3*ls, 5*ls)

%clear; clc; close all;

%% 1. Geometry Setup
% Use spherical frame for infinite medium simulation
geometry = struct( 'dimension', 3 , ...   
                   'frame', 'spherical' ); 
geometry.bnd(1) = struct('dir',3,'val',10);
%% 2. Source Properties
% "Uniform" direction is crucial to observe the diffusion peak (hump)
% at larger distances (e.g., L = 5*ls), matching the paper's isotropic source.
source = struct( 'numberParticles', 5e5 , ...  % High particle count for smooth curves
                 'type', 'point', ...         
                 'position', [0 0 0], ...     
                 'direction', 'uniform', ...  % Isotropic emission (light bulb vs laser)
                 'radial', 0, ...             
                 'lambda', 0.5e-3 );          % Small non-zero width to avoid singularities

%% 3. Observation Configuration
% We need to observe far enough to capture L = 5*ls. 
% If ls ~ 10mm, we need observation up to at least 50mm.
obs_radius_bins = (0.5 : 1.0 : 60.5) * 1e-3; % 1mm spatial resolution

% Time needs to be long enough for the wave to travel 50mm and decay.
% Ballistic time for 50mm is ~110us. We simulate up to 2500us to see the tail.
obs_time = 0 : 2e-6 : 2500e-6; 

% Integration angles for spherical coordinates
% Z is elevation [-pi/2, pi/2], Y is azimuth [-pi, pi]
observation = struct('x', obs_radius_bins, ...      
                     'y', [-pi pi], ...           
                     'z', [-pi/2 pi/2], ...       
                     'directions', [0 pi], ...    
                     'time', obs_time );          

%% 4. Material Properties (Weak Scattering Regime)
freq = 1e6;       
V_paper = 450;    
variance = 0.39;       % From your calibration 
corr_length = 0.0003;  % From your calibration 

material = MaterialClass( geometry, ...
                          freq, ...
                          true, ...          % acoustics = true
                          V_paper, ...             
                          [variance variance], ...     
                          -0.5, ...          
                          'exp', ...         
                          corr_length);            

% Calculate scattering cross-sections (Essential step)
fprintf('Preparing scattering parameters...\n');
material = prepareSigma(material, geometry.dimension);
ls = material.meanFreePath;
fprintf('Calculated Mean Free Path l_s = %.2f mm\n', ls*1000);

%% 5. Run Monte Carlo Simulation
fprintf('Running Simulation (Infinite Medium)... Please wait.\n');
obs = radiativeTransferUnbounded( geometry, source, material, observation );

%% 6. Plotting Results (Figure 1 Reproduction)
% Goal: Plot Energy Density vs Normalized Time for L = 1ls, 3ls, 5ls.
% The curves should separate and peak at t/t_scale = 1, 3, and 5.

ratios = [1, 3, 5]; 
colors = {'k', 'b', 'r'}; % Black, Blue, Red
lines  = {'-', '-', '-'};

% Setup Figure
figure('Name', 'Reproduction of Figure 1 (Zhou et al. 2021)', 'Color', 'w');
axes('YScale', 'log'); % Logarithmic Y-axis is mandatory
hold on; box on; grid on;

% *** Time Normalization Scale ***
% To match the paper, time is normalized by the time required to travel 
% ONE mean free path (ls), not the specific arrival time of each curve.
t_scale = ls / V_paper; 

fprintf('\n--- Plotting Data ---\n');

for i = 1:length(ratios)
    % 1. Determine physical target distance L
    target_L = ratios(i) * ls;
    
    % 2. Find the closest spatial bin index
    bin_centers = (obs_radius_bins(1:end-1) + obs_radius_bins(2:end)) / 2;
    [min_diff, idx] = min(abs(bin_centers - target_L));
    
    % Validation
    if min_diff > 2e-3
        warning('Target distance %.1f mm is far from bin center.', target_L*1000);
    end
    fprintf('Curve %d: Target L = %.1f mm (Ratio %d)\n', i, target_L*1000, ratios(i));

    % 3. Extract Energy Density
    % Squeeze removes singleton dimensions (angles)
    E_curve = squeeze(obs.energyDensity(idx, :, :));
    E_curve = E_curve(:); % Ensure column vector
    
    % 4. Normalize Time Axis
    % t_norm = t / (ls / V)
    t_norm = obs_time(:) / t_scale;
    
    % 5. Plot
    % Match vector lengths to avoid errors
    len = min(length(E_curve), length(t_norm));
    
    semilogy(t_norm(1:len), E_curve(1:len), ...
             'Color', colors{i}, ...
             'LineWidth', 2, ...
             'LineStyle', lines{i}, ...
             'DisplayName', sprintf('L = %d l_s', ratios(i)));
end

% Formatting
xlabel('Normalized Time (t / (l_s/V))');
ylabel('Energy Density (J/m^3)');
title(['Infinite Medium Energy Profiles (l_s \approx ' num2str(ls*1000, '%.1f') ' mm)']);
legend('show', 'Location', 'northeast');

% Set limits to match the paper's aesthetic
xlim([0 25]);       % View range
ylim([1e0 1e7]);    % Energy dynamic range

% Add reference lines for arrival times
xline(1, 'k:', 'HandleVisibility', 'off');
xline(3, 'b:', 'HandleVisibility', 'off');
xline(5, 'r:', 'HandleVisibility', 'off');

hold off;
fprintf('Plot complete.\n');