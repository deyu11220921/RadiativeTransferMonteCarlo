clear all 
%% 1. Geometry (Finite Cylinder)
% This defines a finite-sized cylindrical body
L = 10;          % Cylinder length (thickness)
w = L/4;          % Cylinder diameter (based on paper w = 2.5 mm)
R_cyl = w/2;      % Cylinder radius

geometry = struct( 'dimension', 3 , ...   
                   'frame', 'cylindrical' ); 
% Under 'cylindrical' coordinate system:
% x -> r (radius)
% y -> phi (azimuthal angle)
% z -> z (height)

% Define Boundaries - Assuming Perfect Reflection
% dir=4 represents the cylinder side wall radius (r), dir=3 represents the Z-axis planes
geometry.bnd(1) = struct('dir', 4, 'val', R_cyl); % Side wall at r = 1.25
geometry.bnd(2) = struct('dir', 3, 'val', 0);     % Bottom surface at z = 0
geometry.bnd(3) = struct('dir', 3, 'val', L);     % Top surface at z = 10 (Exit)

%% 2. Point Source
% Source is placed at the center of the bottom surface, emitting inwards
source = struct( 'numberParticles', 1e5 , ... 
                 'type', 'point', ...          
                 'position', [0 0 0], ...     % (x,y,z) in Cartesian
                 'direction', [0 0 1], ...  
                 'lambda', eps ); 

%% 3. Material Properties
% To reproduce Figure 4, we need to set the scattering mean free path (ls)
factor = 10;  % You can modify this to 1, 3, or 10 corresponding to Fig 4 (a)(b)(c)
ls = L / factor; 

freq = 1e6;           
material.v = 450000; % Velocity in mm/s
material.acoustics = true; 

% Scattering coefficient sigma = 1/ls
% Isotropic scattering: 1/4/pi
material.sigma{1} = @(th) 1/4/pi*ones(size(th))*material.v*factor/L; 

% Absorption: Qi = 500
Qi = 500;
omega = 2 * pi * freq;

 material.rho = 2300;
material = prepareSigma(material, geometry.dimension);
fprintf('Setup: L=%.1f, ls=%.2f (L/ls=%.1f)\n', L, ls, factor);

%% 4. Observation (Key modifications)
tb = L / material.v;
dt = material.meanFreeTime/10;
obs_time = 0:dt:20*tb;

% We need to observe the energy at z = L.
% In cylindrical coordinates, we set observation bins for z.
% To locate dz=10 (z=L), we set a grid in the z direction ensuring the last grid includes L.

dz_step = 0.25; % Z-axis observation resolution
z_bins = 0:dz_step:L; % Z-axis grid [0, 0.5, ..., 9.5, 10.0]

% x (radius r): For "On-axis" detection, take a small range near 0.
% For transmission across the entire surface, take 0 to R_cyl.
% Here, a detailed r-grid is set for flexible selection later.
r_bins = 0:0.25:R_cyl; 

observation = struct('x', [0 0.25], ...       % r (Radius)
                     'y', [-pi pi], ...     % Azimuth (integration over all angles)
                     'z', z_bins, ...       % Z (Height) - includes z=10
                     'directions', [0 pi], ...    
                     'time', obs_time );     

%% 5. Solution
fprintf('Running Simulation (Finite Cylinder)...\n');
tic; 
obs = radiativeTransferUnbounded( geometry, source, material, observation ); 
toc;

%% 6. Data Extraction for [20 1 201]
en = squeeze(obs.energyDensity);
y_sim_mc = en(end, :) * 1e9; 

% 2. Process the time axis
% Ensure x_sim length matches y_sim_mc
x_sim = obs.t / tb; 
att = exp(- (omega * observation.time) / Qi);
y_sim_mc = y_sim_mc .* att;

ref1 = Datacylindre(1); % factor=1
ref3 = Datacylindre(2); % factor=3
ref10 = Datacylindre(3); % factor=10

%% 7. Plotting
figure;
h1 = semilogy(x_sim, y_sim_mc, 'k-', 'LineWidth', 2); 
hold on; 
h2 = semilogy(ref10(:,1), ref10(:,2), 'r--', 'DisplayName', 'Reference 10 (Long)');
grid on;

% Add legend and labels
title(['Finite Cylinder (L/l_s = ' num2str(factor) ') - Transmission at z=L']);
xlabel('Normalized Time (t / t_b)');
ylabel('Normalized Energy Density');

% Mark the ballistic arrival time (t/tb = 1)
xline(1, 'k--', 'Ballistic Arrival'); 

% Set axis limits (referencing paper Figure 4)
xlim([0 20]);