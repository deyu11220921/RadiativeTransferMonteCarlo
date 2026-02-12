
clear all 
%% 1. Geometry
% Infinite medium (spherical frame)
geometry = struct( 'dimension', 3 , ...   
                   'frame', 'spherical' ); 
               
%% 2. Point Source
% Small lambda to avoid numerical singularity at origin
source = struct( 'numberParticles', 1e5 , ... 
                 'type', 'point', ...         
                 'position', [0 0 0], ...     
                 'direction', 'uniform',...%[0 0 1], ...    
                 'lambda', eps ); 

%% 3. factor
factor = 3;
L = 10; dr = 0.25; m = ceil(L/dr);
% Spatial bins (Radius R)
if rem(m,2)==0
    obs_radius_bins = L-((m-1):-2:-1)*dr;
else
    obs_radius_bins = L-((m-2):-2:-1)*dr;
end
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
% Fix: Z is elevation in [-pi/2, pi/2], not colatitude [0, pi]

%% 5. observation
tb = L/material.v;
dt = material.meanFreeTime/10;
obs_time = 0:dt:20*tb; 
observation = struct('x', obs_radius_bins, ...      
                     'y', [-pi pi], ...           % Azimuth (Full circle)
                     'z', [-pi/2 pi/2], ...       % Elevation (Full sphere)
                     'directions', [0 pi], ...    
                     'time', obs_time );     


            
%% 6. Solution
fprintf('Running Simulation (Infinite Medium)...\n');
tic ; obs = radiativeTransferUnbounded( geometry, source, material, observation ); toc;

%% 7. Plotting Results
en = squeeze(obs.energyDensity);
y_sim_mc= en(end, :) * 1e9; % 
x_sim = obs.t / tb;       % 

Qi = 500;           
           
omega = 2 * pi * freq ; 
att = exp(- (omega * observation.time) / Qi);
y_sim_mc = y_sim_mc .* att;

% --- 1. Load specific datasets ---
ref1 = Datainifinie(1); %  factor=1
ref3 = Datainifinie(2); %  factor=3
ref10 = Datainifinie(3); % factor=10

h1 = semilogy(x_sim, y_sim_mc, 'k-', 'LineWidth', 2); 
hold on; 
h2 = semilogy(ref3(:,1), ref3(:,2), 'r--', 'DisplayName', 'Reference 1 (Long)');

% --- 3. Figure Decorations ---
grid on;
xlabel('Normalized Time (t / t_b)');
ylabel('Energy Density (J/m^3)');
title(['Infinite Medium']);
% Add the ballistic arrival reference line
xline(1, 'r--', 'Ballistic Arrival', 'LabelVerticalAlignment', 'bottom');
% Set axes limits
xlim([0 20]);
% Add a legend to distinguish the two lines
legend([h1, h2], {'Monte Carlo', 'Reference Data'}, 'Location', 'northeast');
hold off; 