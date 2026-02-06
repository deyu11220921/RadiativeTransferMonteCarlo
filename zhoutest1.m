

%% 1. Geometry
% Infinite medium (spherical frame)
geometry = struct( 'dimension', 3 , ...   
                   'frame', 'spherical' ); 
geometry.bnd(1) = struct('dir',3,'val',10);
%% 2. Point Source
% Small lambda to avoid numerical singularity at origin
source = struct( 'numberParticles', 1e3 , ... 
                 'type', 'point', ...         
                 'position', [0 0 0], ...     
                 'direction', 'uniform',...%[0 0 1], ...    
                 'radial', 0, ...             
                 'lambda', 1e-4 ); 

%% 3. Observation
% Spatial bins (Radius R)
obs_radius_bins = (0.5 : 1.0 : 20.5) * 1e-3; 

% Time vector (Simulate up to 500us)
obs_time = 0 : 0.5e-6 : 500e-6; 

% Fix: Z is elevation in [-pi/2, pi/2], not colatitude [0, pi]
observation = struct('x', obs_radius_bins, ...      
                     'y', [-pi pi], ...           % Azimuth (Full circle)
                     'z', [-pi/2 pi/2], ...       % Elevation (Full sphere)
                     'directions', [0 pi], ...    
                     'time', obs_time );          

%% 4. Material Properties
freq = 1e6;       
V_paper = 450;    
variance = 0.39;        % From calibration   1(0.123)   3(0.2153)  10(0.39) 
corr_length = 0.0003;  % From calibration    0.0003

material = MaterialClass( geometry, ...
                          freq, ...
                          true, ...          % acoustics = true
                          V_paper, ...             
                          [variance variance], ...     
                          -0.5, ...          
                          'exp', ...         
                          corr_length);            

% Prepare scattering matrices
fprintf('Preparing scattering cross-sections...\n');
material = prepareSigma(material, geometry.dimension);
fprintf('Mean Free Path l_s = %.2f mm\n', material.meanFreePath*1000);

%% 5. Solution
fprintf('Running Simulation (Infinite Medium)...\n');
obs = radiativeTransferUnbounded( geometry, source, material, observation );

%% 6. Plotting Results
target_L = 10e-3; % Target distance: 10 mm

% Locate the spatial bin for R = 10mm
bin_centers = (obs_radius_bins(1:end-1) + obs_radius_bins(2:end)) / 2;
[~, idx_L] = min(abs(bin_centers - target_L)); 

% Extract Energy Density
% obs.energyDensity dimensions: [Nx, 1, Nt] -> Squeeze to [Nx, Nt]
E_total = squeeze(obs.energyDensity(idx_L, :, :));

% Normalize axes
tb = target_L / V_paper; % Ballistic time
t_norm = obs_time / tb;

% Plot (Semilogy for energy decay)
figure;
semilogy(t_norm, E_total, 'k-', 'LineWidth', 2);
grid on;
xlabel('Normalized Time (t / t_b)');
ylabel('Energy Density (J/m^3)');
title(['Infinite Medium (L = 10mm, l_s \approx 10mm)']);
xline(1, 'r--', 'Ballistic Arrival');
xlim([0 20]);