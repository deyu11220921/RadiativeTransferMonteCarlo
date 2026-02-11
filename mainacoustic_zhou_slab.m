%% Zhou 2021 Figure 4(b) Reproduction - Finite Slab (Thin Plate) Case
% Target: Anomalous Transport (Ballistic & Diffusion Peaks Merge)

clear; clc; close all;

%% 1. Geometry (Slab = Wide Cylinder)

L_sample = 10e-3;        
r_sample = 20e-3; % r = W/2 = 20mm

geometry = struct('dimension', 3, ...
                  'frame', 'cylindrical'); 

% (Perfect Boundaries) (Reverberations)
% dir=3 -> Z, dir=4 -> R
geometry.bnd(1) = struct('dir', 3, 'val', 0);          % Bottom (Z=0) 
geometry.bnd(2) = struct('dir', 3, 'val', L_sample);   % Top (Z=L) 
geometry.bnd(3) = struct('dir', 4, 'val', r_sample);   % Lateral (R=20mm)

%% 2. Point Source

source = struct('numberParticles', 1e5, ... % 
                'type', 'point', ...
                'position', [0 0 0], ...    % 
                'direction', [0 0 1], ...   % 
                'radial', 0, ...
                'lambda', 1e-4); 

%% 3. Observation (Dimension Correction - Bypassing Library Bug) 
% At least one dimension must have more than 1 bin.
% Otherwise, 'initializeObservation' will fail with "Matrix dimensions must agree" because 'ibins' is empty.
% 
% Solution: We forcibly split the Y (Phi) axis into two bins ([-pi, 0] and [0, pi]).
% This tricks the library into recognizing a valid 2D observation, allowing it to run normally.
% We will simply sum/average these two symmetric parts during post-processing.

% Z-axis (Depth): Observe transmitted energy at Z=L
dz = 0.1e-3; 
obs_z_bins = [L_sample-dz, L_sample+dz]; 

% X-axis (Radius): Only observe on-axis region [0, 2mm]
detector_radius = 2e-3; 
obs_r_bins = [0 detector_radius]; 

% [Key Correction] Y-axis (Phi): 
% The previous setting [-Inf Inf] caused two issues: 1. Crash, 2. Energy density divided by Inf becomes 0.
% Here we change it to [-pi 0 pi] to create 2 bins, activating the library's histogram logic.
obs_phi = [-pi 0 pi]; 

% Time: Normalized time t/tb
V_wave = 450; 
tb = L_sample / V_wave; 
obs_time = 0 : (tb/50) : (20 * tb); 

observation = struct('x', obs_r_bins, ...    % Dimension 1: Radius
                     'y', obs_phi, ...       % Dimension 2: Phi (Dummy split for stability)
                     'z', obs_z_bins, ...    % Dimension 3: Depth
                     'directions', [0 pi], ... 
                     'time', obs_time);
%% 4. Material Properties 
freq = 1e6;        
acoustics = true; 
corr_func = 'exp'; 

% =======================================================


calibrated_variance = 0.39;   % (L/ls = 1, 3, or 10) 1(0.123)   3(0.2153)  10(0.39)
calibrated_corr_len = 0.0003;  % 0.0003
% =======================================================

material = MaterialClass(geometry, freq, acoustics, V_wave, ...
                         [calibrated_variance calibrated_variance], ...
                         -0.5, corr_func, calibrated_corr_len);

material.rho = 2300;
fprintf('Preparing Scattering Matrix...\n');
material = prepareSigma(material, geometry.dimension);
ls_actual = material.meanFreePath;

fprintf('\n---------------- CHECK PARAMETERS ----------------\n');
fprintf('Target Ratio L/l_s = 3.0\n');
fprintf('Calculated l_s     = %.4f mm\n', ls_actual*1000);
fprintf('Actual Ratio L/l_s = %.2f\n', L_sample / ls_actual);
fprintf('--------------------------------------------------\n\n');

%% 5. Simulation (Running Code)
fprintf('Running Simulation (Finite Slab - On Axis Detection)...\n');
obs = radiativeTransferUnbounded(geometry, source, material, observation);

%% 6. Post-Processing & Plotting (Corrected - FIXED DIMENSION BUG)
% This section processes the raw simulation data to match Zhou 2021 Figure 4(b).

% 1. Extract Raw Data
if isfield(obs, 't')
    time_vec = obs.t;
else
    time_vec = obs.time;
end
time_vec = time_vec(:)'; % Ensure row vector

% --- [CRITICAL FIX START] ---
% obs.energyDensity dimensions: [Radius(1) x Phi(2) x Time(Nt)]
% We must ONLY sum over spatial dimensions (1 and 2).
% DO NOT sum over dimension 3 (Time), or you will destroy the signal!

E_sum_spatial = sum(sum(obs.energyDensity, 1), 2); 
% Result is now [1 x 1 x Nt]

E_raw = squeeze(E_sum_spatial); 
% Result is now [Nt x 1] (or [1 x Nt] depending on Matlab version)
E_raw = E_raw(:)'; % Force to row vector [1 x Nt] matches time_vec
% --- [CRITICAL FIX END] ---

% 2. Introduce Intrinsic Absorption Qi (Qi=500)
Qi = 500;           
f = 1e6;            
omega = 2 * pi * f; 

% Calculate attenuation factor: exp(-omega * t / Qi)
absorption_factor = exp(- (omega * time_vec) / Qi);

% Apply absorption 
E_final = E_raw .* absorption_factor; 

% 3. Normalize Time Axis
V_wave = 450;       
L_sample = 10e-3;   
tb = L_sample / V_wave; 
t_norm = time_vec / tb;

% 4. Smoothing (DISABLED)
E_plot = E_final; 

% 5. Plotting (Target: Figure 4b)
figure('Color', 'w', 'Name', 'Zhou Fig 4b Reproduction (Slab)');

% Plot
semilogy(t_norm, E_plot, 'b.-', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', 'Slab On-axis (Simulated)');

grid on;
xlabel('Normalized Time (t / t_b)', 'FontSize', 12);
ylabel('Energy Density (J/m^3)', 'FontSize', 12);

% Add reference line for Ballistic Arrival
xline(1, 'r--', 'Ballistic Arrival', 'LabelVerticalAlignment', 'bottom');

% Adjust Axes
xlim([0 20]); 

% Smart Y-Limits
max_E = max(E_plot);
if max_E > 0
    ylim([max_E * 1e-4, max_E * 2]); 
else
    warning('Signal is zero. Check source/detector alignment.');
end

title({sprintf('Finite Slab On-Axis (Q_i=%d)', Qi); ...
       'Expectation: Sharp Ballistic Peak at t=1 + Reverberations'}, 'FontSize', 12);

legend('show');