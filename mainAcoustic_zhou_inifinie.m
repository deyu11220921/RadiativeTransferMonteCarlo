

%% 1. Geometry
% Infinite medium (spherical frame)
geometry = struct( 'dimension', 3 , ...   
                   'frame', 'spherical' ); 
%% 2. Point Source
% Small lambda to avoid numerical singularity at origin
source = struct( 'numberParticles', 2e5 , ... 
                 'type', 'point', ...         
                 'position', [0 0 0], ...     
                 'direction', 'uniform',...%[0 0 1], ...    
                 'lambda', eps ); 

%% 3. Observation
factor = 10;
L = 10; dr = 0.25; m = ceil(L/dr);
% Spatial bins (Radius R)
if rem(m,2)==0
    obs_radius_bins = L-((m-1):-2:-1)*dr;
else
    obs_radius_bins = L-((m-2):-2:-1)*dr;
end
tb = L/material.v;
dt = material.meanFreeTime;
obs_time = 0:dt:20*tb; 

% Fix: Z is elevation in [-pi/2, pi/2], not colatitude [0, pi]
observation = struct('x', obs_radius_bins, ...      
                     'y', [-pi pi], ...           % Azimuth (Full circle)
                     'z', [-pi/2 pi/2], ...       % Elevation (Full sphere)
                     'directions', [0 pi], ...    
                     'time', obs_time );          

tic; obs = radiativeTransferUnbounded( geometry, source, material, observation ); toc;
                 
%% 4. Material Properties
freq = 1e6;          
material.v = 450000; 
material.acoustics = true; 
material.sigma{1} = @(th) 1/4/pi*ones(size(th))*material.v*factor/L;

material = prepareSigma(material,geometry.dimension);
% Prepare scattering matrices
fprintf('Preparing scattering cross-sections...\n');
material = prepareSigma(material, geometry.dimension);
fprintf('Mean Free Path l_s = %.2f mm\n', material.meanFreePath);

%% 5. Solution
fprintf('Running Simulation (Infinite Medium)...\n');
obs = radiativeTransferUnbounded( geometry, source, material, observation );

%% 6. Plotting Results
en = squeeze(obs.energyDensity);
y_sim = en(end, :) * 1e9; % 
x_sim = obs.t / tb;       % 

Qi = 500;           
f = 1e6;             
omega = 2 * pi * f; 
y_sim = y_sim .* exp(- (omega * obs.t) / Qi); 

figure;
semilogy(x_sim, y_sim, 'k-', 'LineWidth', 2);
grid on;
xlabel('Normalized Time (t / t_b)');
ylabel('Energy Density (J/m^3)');
title(['Infinite Medium (L = 10mm, l_s \approx 10mm)']);
xline(1, 'r--', 'Ballistic Arrival');
xlim([0 20]);