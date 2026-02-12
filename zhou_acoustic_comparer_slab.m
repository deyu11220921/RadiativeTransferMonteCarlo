clear all 
%% 1. Geometry (Infinite Slab / Thin Plate)
% Paper parameters: Slab geometry, L = 10 mm. 
% Laterally infinite; boundaries only exist at the top and bottom surfaces.
L = 10;          
w=4*L;
R_cyl = w/2;      % Cylinder radius

geometry = struct( 'dimension', 3 , ...   
                   'frame', 'cylindrical' ); 
% Using cylindrical coordinates (r, phi, z) to facilitate on-axis observation
% x -> r (radius)
% y -> phi (azimuthal angle)
% z -> z (depth/height)

% --- Boundary Definitions ---
% A Slab only has boundaries in the Z direction (z=0 and z=L).
% No lateral side boundaries (dir=4 removed).
geometry.bnd(1) = struct('dir', 4, 'val', R_cyl); % Side wall at r = 1.25
geometry.bnd(2) = struct('dir', 3, 'val', 0);     % Bottom surface at z = 0
geometry.bnd(3) = struct('dir', 3, 'val', L);     % Top surface at z = 10 (Exit)

%% 2. Point Source
% Increased particle count to obtain a smoother on-axis signal
source = struct( 'numberParticles', 2e6 , ... 
                 'type', 'point', ...          
                 'position', [0 0 0], ...     
                 'direction', [0 0 1], ...  % Paper Fig.2 describes point source as isotropic
                 'lambda', eps ); 

%% 3. Material Properties
% Set the scattering mean free path (ls)
% Paper Figure 4 compares factor = 1, 3, 10
factor = 1;  % Modify this to 1, 3, or 10 to reproduce different curves
ls = L / factor; 

freq = 1e6;           
material.v = 450000; % Velocity in mm/s
material.acoustics = true; 

% Define scattering cross-section sigma = 1/ls
% Isotropic scattering: 1/4/pi
material.sigma{1} = @(th) 1/4/pi*ones(size(th))*material.v*factor/L; 

% Absorption
Qi = 500;
omega = 2 * pi * freq;

material = prepareSigma(material, geometry.dimension);
fprintf('Setup: Slab Geometry, L=%.1f, ls=%.2f (L/ls=%.1f)\n', L, ls, factor);

%% 4. Observation (On-axis Detection at z=L)
tb = L / material.v;
dt = material.meanFreeTime/10; 
obs_time = 0:dt:20*tb; 

% --- Key Modification: Observation Regions ---
% 1. Z-direction: We need the energy exactly at the exit point z = L.
%    Take a thin layer in [L-dz, L] as the detection region.
dz_step = 0.25; 
z_bins = 0:dz_step:L;

% 2. R-direction (x): To reproduce the "On-axis" effect, r must be restricted.
r_detect_radius = 2.0; 
r_bins = [0 r_detect_radius]; 

observation = struct('x', [0 2.5], ...         % r (On-axis detector)
                     'y', [-pi pi], ...       % Azimuth (integration over all angles)
                     'z', z_bins, ...         % Z (at exit surface)
                     'directions', [0 pi], ...    
                     'time', obs_time );     

%% 5. Solution
fprintf('Running Simulation (Infinite Slab)...\n');
tic; 
obs = radiativeTransferUnbounded( geometry, source, material, observation ); 
toc;

%% 6. Data Extraction
en = squeeze(obs.energyDensity);
y_sim_mc = en(end, :) * 1e9; 

% 2. Process the time axis
% Ensure x_sim length matches y_sim_mc
x_sim = obs.t / tb; 
att = exp(- (omega * observation.time) / Qi);
y_sim_mc = y_sim_mc .* att;

% Reference data loaders
ref1 = Dataslab(1); % factor=1
ref3 = Dataslab(2); % factor=3
ref10 = Dataslab(3); % factor=10

%% 7. Plotting
figure;
h1 = semilogy(x_sim, y_sim_mc, 'b-', 'LineWidth', 2, 'DisplayName', ['Slab MC (L/ls=' num2str(factor) ')']); 
hold on; 
h2 = semilogy(ref1(:,1), ref1(:,2), 'r--', 'DisplayName', 'Reference 10 (Long)');

% Add Ballistic Arrival reference line
xline(1, 'k--', 'Ballistic Arrival'); 

grid on;
title(['Slab Geometry (On-axis detection) - L/l_s = ' num2str(factor)]);
xlabel('Normalized Time (t / t_b)');
ylabel('Energy Density (J/m^3)');
xlim([0 20]);
legend('Location', 'northeast');
hold off;